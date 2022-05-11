/* flagstats-tbb.cc -- parallel bam flagstats using mmbam and tbb
 *
 * Author: Yi Qiao <yi.qiao@genetics.utah.edu>
 *
 * This code example is adopted from the bam_stats.c in the samtools project.
 * It demonstrates how to use libmmbam to compute sequence alignment stats
 * in parallel using the Intel Thread Building Block library for multi-
 * threading. 
 *
 * The original copy right claims, author information, and LICENSE is included
 * below
 */

/*  bam_stat.c -- flagstat subcommand.

    Copyright (C) 2009, 2011, 2013-2015, 2019, 2021 Genome Research Ltd.

    Author: Heng Li <lh3@sanger.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <fstream>
#include <string>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <mmbam/mfile.h>
#include <mmbam/bam.h>
#include <mmbam/index.h>

// The following defs are taken from htslib/sam.h
/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256
/*! @abstract QC failure */
#define BAM_FQCFAIL      512
/*! @abstract optical or PCR duplicate */
#define BAM_FDUP        1024
/*! @abstract supplementary alignment */
#define BAM_FSUPPLEMENTARY 2048


typedef struct {
    long long n_reads[2], n_mapped[2], n_pair_all[2], n_pair_map[2], n_pair_good[2];
    long long n_sgltn[2], n_read1[2], n_read2[2];
    long long n_dup[2];
    long long n_diffchr[2], n_diffhigh[2];
    long long n_secondary[2], n_supp[2];
    long long n_primary[2], n_pmapped[2], n_pdup[2];
} bam_flagstat_t;

bam_flagstat_t *bam_flagstat_add(bam_flagstat_t* op1, bam_flagstat_t* op2) {
    long long *op1_view = (long long *)op1;
    long long *op2_view = (long long *)op2;
    for(int i=0; i<32; i++) op1_view[i] += op2_view[i];

    return op1;
}

template<class BAM_T>
inline static void flagstat_loop(bam_flagstat_t *s, const BAM_T& b)
{
    auto *c = &b;
    int w = (c->flag & BAM_FQCFAIL)? 1 : 0;
    if(w) throw std::runtime_error("fail qc");
    ++s->n_reads[w];
    if (c->flag & BAM_FSECONDARY ) {
        ++s->n_secondary[w];
    } else if (c->flag & BAM_FSUPPLEMENTARY ) {
        ++s->n_supp[w];
    } else {
        ++s->n_primary[w];

        if (c->flag & BAM_FPAIRED) {
            ++s->n_pair_all[w];
            if ((c->flag & BAM_FPROPER_PAIR) && !(c->flag & BAM_FUNMAP) ) ++s->n_pair_good[w];
            if (c->flag & BAM_FREAD1) ++s->n_read1[w];
            if (c->flag & BAM_FREAD2) ++s->n_read2[w];
            if ((c->flag & BAM_FMUNMAP) && !(c->flag & BAM_FUNMAP)) ++s->n_sgltn[w];
            if (!(c->flag & BAM_FUNMAP) && !(c->flag & BAM_FMUNMAP)) {
                ++s->n_pair_map[w];
                if (c->next_ref_id != c->ref_id) {
                    ++s->n_diffchr[w];
                    if (c->mapq>= 5) ++s->n_diffhigh[w];
                }
            }
        }

        if (!(c->flag & BAM_FUNMAP)) ++s->n_pmapped[w];
        if (c->flag & BAM_FDUP) ++s->n_pdup[w];
    }
    if (!(c->flag & BAM_FUNMAP)) ++s->n_mapped[w];
    if (c->flag & BAM_FDUP) ++s->n_dup[w];
}

template<class MFILE_T, class REGIONS_T>
bam_flagstat_t *bam_flagstat_core(MFILE_T& mfile, const REGIONS_T& regions)
{
    bam_flagstat_t *s;

    s = (bam_flagstat_t*)calloc(regions.size(), sizeof(bam_flagstat_t));


    parallel_for(tbb::blocked_range<long>(0, regions.size()), [&](tbb::blocked_range<long>& r) {
        for(long i=r.begin(); i<r.end(); ++i) {
            auto bam_records = bam_load_block(mfile, regions[i].first, regions[i].second);
            auto bam_it_beg = bam_iterator(reinterpret_cast<const bam_rec_t*>(&bam_records[0]));
            auto bam_it_end = bam_iterator(reinterpret_cast<const bam_rec_t*>(&bam_records[0] + bam_records.size()));
            auto bam_it = bam_it_beg;
            while(bam_it < bam_it_end) {
                flagstat_loop(s+i, *bam_it);
                bam_it++;
            }
        }
    });

    bam_flagstat_t *ts = (bam_flagstat_t*)calloc(1, sizeof(bam_flagstat_t));
    memset(ts, 0, sizeof(bam_flagstat_t));
    for(long i=0; i<regions.size(); i++) ts = bam_flagstat_add(ts, s+i);

    free(s);

    return ts;
}

static const char *percent(char *buffer, long long n, long long total)
{
    if (total != 0) sprintf(buffer, "%.2f%%", (float)n / total * 100.0);
    else strcpy(buffer, "N/A");
    return buffer;
}

static void usage_exit(FILE *fp, int exit_status)
{
    fprintf(fp, "Usage: samtools flagstat [options] <in.bam>\n");
    fprintf(fp, "  -O, --");
    fprintf(fp, "output-fmt FORMAT[,OPT[=VAL]]...\n"
            "               Specify output format (json, tsv)\n");
    exit(exit_status);
}

static void out_fmt_default(bam_flagstat_t *s)
{
    char b0[16], b1[16];
    printf("%lld + %lld in total (QC-passed reads + QC-failed reads)\n", s->n_reads[0], s->n_reads[1]);
    printf("%lld + %lld primary\n", s->n_primary[0], s->n_primary[1]);
    printf("%lld + %lld secondary\n", s->n_secondary[0], s->n_secondary[1]);
    printf("%lld + %lld supplementary\n", s->n_supp[0], s->n_supp[1]);
    printf("%lld + %lld duplicates\n", s->n_dup[0], s->n_dup[1]);
    printf("%lld + %lld primary duplicates\n", s->n_pdup[0], s->n_pdup[1]);
    printf("%lld + %lld mapped (%s : %s)\n", s->n_mapped[0], s->n_mapped[1], percent(b0, s->n_mapped[0], s->n_reads[0]), percent(b1, s->n_mapped[1], s->n_reads[1]));
    printf("%lld + %lld primary mapped (%s : %s)\n", s->n_pmapped[0], s->n_pmapped[1], percent(b0, s->n_pmapped[0], s->n_primary[0]), percent(b1, s->n_pmapped[1], s->n_primary[1]));
    printf("%lld + %lld paired in sequencing\n", s->n_pair_all[0], s->n_pair_all[1]);
    printf("%lld + %lld read1\n", s->n_read1[0], s->n_read1[1]);
    printf("%lld + %lld read2\n", s->n_read2[0], s->n_read2[1]);
    printf("%lld + %lld properly paired (%s : %s)\n", s->n_pair_good[0], s->n_pair_good[1], percent(b0, s->n_pair_good[0], s->n_pair_all[0]), percent(b1, s->n_pair_good[1], s->n_pair_all[1]));
    printf("%lld + %lld with itself and mate mapped\n", s->n_pair_map[0], s->n_pair_map[1]);
    printf("%lld + %lld singletons (%s : %s)\n", s->n_sgltn[0], s->n_sgltn[1], percent(b0, s->n_sgltn[0], s->n_pair_all[0]), percent(b1, s->n_sgltn[1], s->n_pair_all[1]));
    printf("%lld + %lld with mate mapped to a different chr\n", s->n_diffchr[0], s->n_diffchr[1]);
    printf("%lld + %lld with mate mapped to a different chr (mapQ>=5)\n", s->n_diffhigh[0], s->n_diffhigh[1]);
}

/*
 * Select flagstats output format to print.
 */
static void output_fmt(bam_flagstat_t *s, const char *out_fmt)
{
 
    out_fmt_default(s);
}

int main(int argc, char *argv[])
{
    bam_flagstat_t *s;
    const char *out_fmt = "default";
    int c, status = EXIT_SUCCESS;

    enum {
        INPUT_FMT_OPTION = CHAR_MAX+1,
    };

    int n = tbb::task_scheduler_init::default_num_threads();
    if(argc>2) n = atoi(argv[2]);

    tbb::task_scheduler_init scheduler(n);
    
    auto mfile = mfile_open(argv[1]);
    auto index = index_read(std::ifstream(std::string(argv[1]) + ".bai"));

    /*header = sam_hdr_read(fp);
    if (header == NULL) {
        fprintf(stderr, "Failed to read header for \"%s\"\n", argv[optind]);
        return 1;
    }*/

    auto regions = index_to_regions(index, mfile->size);

    index_free(index);


    //s = bam_flagstat_core(fp, header);
    s = bam_flagstat_core(mfile, regions);
    if (s) {
        output_fmt(s, out_fmt);
        free(s);
    }

    //sam_hdr_destroy(header);
    //sam_close(fp);
    //sam_global_args_free(&ga);
    return status;
}
