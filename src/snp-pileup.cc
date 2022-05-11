/* snp-pileup.cc -- parallel snp-pileup using mmbam and tbb
 *
 * Author: Yi Qiao <yi.qiao@genetics.utah.edu>
 *
 * This code example reimplements the snp-pileup from the copy number variation
 * calling algorithm FACETS in a functional equivalent manner. It demonstrates
 * how to perform bam file processing in a parallelized manner over the pileup
 * of many locations.
 *
 * Some code snippets are taken from the original snp-pileup.cpp. Since there 
 * is no copyright announcements in the original repository, we wish to declare
 * that the authorship and credits of these snippets go to the orignal authors. 
 */

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <list>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <map>

#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>

#include <mmbam/mbgzf.h>
#include <mmbam/index.h>
#include <mmbam/bam.h>
#include <mmbam/mpileup.h>

template<typename CLK> void timestamp(CLK& clk_start, const std::string& label) {
    auto clk_now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(clk_now - clk_start).count();
    auto secs = duration / 1000;
    auto mils = duration % 1000;
    std::cerr << "["<<secs<<"."<<mils<<"] " << label << std::endl;
}

struct vcf_record {
    std::string chrom;
    int32_t pos;
    uint8_t ref;
    uint8_t alt;

    int ref_count[2] = {0, 0};
    int alt_count[2] = {0, 0};
    int del_count[2] = {0, 0};
    int err_count[2] = {0, 0};

    int depth[2];
};

void load_vcf(std::vector<vcf_record>& records, std::vector<std::pair<size_t, size_t>>& work_ranges, const char *filename, const indices_t& indices, const std::vector<std::map<std::string, uint32_t>>& tid_maps);

constexpr uint16_t flag_fail = READ_FAILED_QUALITY_CHECKS | READ_SECONDARY_ALIGNMENT | READ_PCR_OPTICAL_DUPLICATE | READ_UNMAPPED;

bool filter_predicate(const bam_rec_t& bam_rec) {
    if (bam_rec.mapq < 1) return false;
    if (bam_rec.flag & flag_fail) return false;
    if (bam_rec.flag & READ_PAIRED && !(bam_rec.flag & READ_MAPPED_PROPER_PAIR)) return false;

    return true;
}

void snp_mpileup(std::vector<vcf_record>& vcf_records,
        const mfiles_t& mfiles, const indices_t& indices,
        const std::vector<std::map<std::string, uint32_t>>& chr_tid_maps,
        bool has_prefix, size_t batch_size = 10000) {

    tbb::parallel_for(tbb::blocked_range<size_t>(0, vcf_records.size(), batch_size), [&](tbb::blocked_range<size_t> r) {

        size_t batch_first_idx = r.begin();
        size_t batch_end_idx = r.end();

        // as long as we haven't exausted our batch
        while(batch_first_idx < batch_end_idx) {
            auto curr_vcf  = vcf_records.begin() + batch_first_idx;
            auto last_vcf  = vcf_records.begin() + batch_end_idx - 1;
            // TODO: following can be replaced with binary search
            while((curr_vcf->chrom != last_vcf->chrom || last_vcf->pos - curr_vcf->pos > 65536) && curr_vcf != last_vcf) {
                batch_end_idx--;
                last_vcf = vcf_records.begin() + (batch_end_idx - 1);
            }

            // from curr_vcf to last_vcf should be on the same chr
            batch_first_idx = batch_end_idx; // setup next iteration

            // reset batch-end
            batch_end_idx = r.end(); 

            auto batch_pos_start = curr_vcf->pos;
            auto batch_pos_end = last_vcf->pos;

            // pileup from curr_vcf to last_vcf
            mpileup(mfiles, indices, chr_tid_maps[0].at(curr_vcf->chrom),
                    curr_vcf->pos, last_vcf->pos + 1, filter_predicate,
                    [&vcf_records, &curr_vcf, &last_vcf, &batch_pos_start, &batch_pos_end](const auto& p) {
                if(curr_vcf > last_vcf) return false;


                // catch up vcf
                while(p.pos > curr_vcf->pos && curr_vcf < last_vcf) curr_vcf++;

                if(p.pos == curr_vcf->pos) {
                    // update curr_vcf
                    for(size_t i_file=0; i_file<2; i_file++) {
                        curr_vcf->depth[i_file] = p.depth[i_file];
                        
                        auto* buffer_start = p.reads_buffer[i_file]->data();

                        for(size_t i=0; i<p.depth[i_file]; i++) {
                            auto& info = p.get_info(i_file, i);
                            const bam_rec_t * bam_rec = BAMREF( buffer_start + info.offset);

                            if(info.is_deletion) curr_vcf->del_count[i_file]++;
                            else {
                                if(bam_bqual_ptr(bam_rec)[info.qpos] < 1) continue;
                                auto seq = bam_seq_ptr(bam_rec);
                                uint8_t base = bam_unpack_base(seq, info.qpos);
                                if(base == curr_vcf->ref) curr_vcf->ref_count[i_file]++;
                                else if(base == curr_vcf->alt) curr_vcf->alt_count[i_file]++;
                                else curr_vcf->err_count[i_file]++;
                            }
                        }
                    }
                    curr_vcf++;
                }

                if(p.ref_id == 20 && p.pos == 10808236 ) {
                    std::cerr<<"visited in region "<<batch_pos_start<<" - "<<batch_pos_end<<std::endl;
                    std::cerr<<"R: "<<curr_vcf->ref_count[0]<<","<<curr_vcf->ref_count[1]<<std::endl;;
                }

            return true;
            });
        }
    });

    
}

std::map<std::string, uint32_t> parse_header(const mfile_t::ptr_t& mfile, const index_t& index, bool& bam_has_prefix) {
    // read header
    bgzf_mfile_proxy_t bgzf_proxy(mfile);
    std::vector<uint8_t> bam_buffer;
    for (auto& bgzf_block : bgzf_proxy) {
        auto de_buff = bgzf_inflate(bgzf_block);
        bam_buffer.insert(bam_buffer.end(), de_buff.cbegin(), de_buff.cend());
        if(bam_buffer_contains_header(bam_buffer)) break;
    }

    const bam_header_t *hdr = 
        reinterpret_cast<const bam_header_t*>(bam_buffer.data());

    uint8_t *buffer_ptr = bam_buffer.data();
    buffer_ptr += sizeof(char) * 4 + sizeof(uint32_t) + sizeof(char) * hdr->l_text;
    uint32_t n_ref = *(uint32_t *)(buffer_ptr);
    buffer_ptr += sizeof(uint32_t); // now bufer_ptr points at first ref

    bam_has_prefix = strncmp((char *)(buffer_ptr + sizeof(uint32_t)), "chr", 3) == 0;

    std::map<std::string, uint32_t> chr_tid_map;

    for (uint32_t i_ref = 0; i_ref < n_ref; i_ref++) {
        uint32_t l_name = *(uint32_t *)(buffer_ptr);
        buffer_ptr += sizeof(uint32_t);
        std::string ref_name(buffer_ptr + (3 * bam_has_prefix), buffer_ptr + l_name - 1);
        buffer_ptr += l_name + sizeof(uint32_t);
        chr_tid_map[ref_name] = i_ref;
    }

    return chr_tid_map;
}

int main(int argc, char** argv) {
    // main logic
    // 1. read all pileup locations from vcf file
    // 1.1. filter out locations with more than 2 alleles
    // 1.2. filter out locations that are not SNPs
    // 2. for each location, pileup each bam file given on the command line
    // 2.1. reads are filtered according to the following criteria 
    // 2.1.1. flag should not indicate fail QC
    // 2.1.2. flag should not indicate secondary alignment
    // 2.1.3. flag should not indicate PCR or optical duplicate
    // 2.1.4. flag should not indicate read unmapped
    // 2.1.5. read tid should not be less than 0
    // 2.1.6. mapping quality should not be less than min_mq (1)
    // 2.1.7. flag should not indicate proper pair if it's paired-end reads
    // 2.2. for each location in each bam file, examine the pileup
    // 2.2.1. ignore locations where total depth is less than min_rc (0)
    // 2.2.1. ignore read if base quality is less than min_bq (1)
    // 2.2.2. count the number of deletions, refs, alts, or others (error)
    // 3. produce an output with the following format
    // 3.1. Chromosome,Position,Ref,Alt,File1R,File1A,File1E,File1D,...
    
    auto clk_start = std::chrono::high_resolution_clock::now();
    timestamp(clk_start, "started");

    // linear pileup
    auto mfile1 = mfile_open(argv[2]);
    auto mfile2 = mfile_open(argv[3]);

    auto index1 = index_read(std::ifstream(std::string(argv[2]) + ".bai"));
    auto index2 = index_read(std::ifstream(std::string(argv[3]) + ".bai"));
    
    mfiles_t mfiles{mfile1, mfile2};
    indices_t indices{index1, index2};
    
    std::vector<std::map<std::string, uint32_t>> chr_tid_maps(2);
    bool bam_has_prefix;
    chr_tid_maps[0] = parse_header(mfile1, index1, bam_has_prefix);
    chr_tid_maps[1] = parse_header(mfile2, index2, bam_has_prefix);

    std::cout<<"Chromosome,Position,Ref,Alt,File1R,File1A,File1E,File1D,File2R,File2A,File2E,File2D";
    std::cout<<std::endl;

    timestamp(clk_start, "header parsed");

    std::vector<vcf_record> records;
    std::vector<std::pair<size_t, size_t>> work_ranges;

    load_vcf(records, work_ranges, argv[1], indices, chr_tid_maps);

    timestamp(clk_start, "vcf loaded, "+std::to_string(records.size()) + " vars, "+std::to_string(work_ranges.size())+" ranges");

    //auto work_ranges = vcf_parallel_ranges(records, indices, chr_tid_maps);
    //timestamp(clk_start, "range calculated, "+std::to_string(work_ranges.size())+" ranges");
    //return 0;

    tbb::parallel_for_each(work_ranges.cbegin(), work_ranges.cend(), [&](auto& r){

        auto curr_vcf = records.begin() + r.first;
        auto last_vcf = records.begin() + (r.second - 1);
        auto end_vcf = records.begin() + r.second;

        // pileup from curr_vcf to last_vcf
        mpileup(mfiles, indices, chr_tid_maps[0].at(curr_vcf->chrom),
                curr_vcf->pos, last_vcf->pos + 1, filter_predicate,
                [&curr_vcf, &end_vcf](const auto& p) {

            // ignore earlier positions
            if (p.pos < curr_vcf->pos) return true;

            // catch up vcf
            while (curr_vcf < end_vcf) {
                if (p.pos > curr_vcf->pos) {
                    curr_vcf++;
                }
                else break;
            }

            if (curr_vcf == end_vcf) return false;

            if (p.pos == curr_vcf->pos) {
                // update curr_vcf
                for (size_t i_file=0; i_file<2; i_file++) {
                    curr_vcf->depth[i_file] = p.depth[i_file];
                    
                    auto* buffer_start = p.reads_buffer[i_file]->data();

                    for (size_t i=0; i<p.depth[i_file]; i++) {
                        auto& info = p.get_info(i_file, i);
                        const bam_rec_t * bam_rec = BAMREF( buffer_start + info.offset);

                        if (info.is_deletion) curr_vcf->del_count[i_file]++;
                        else {
                            if (bam_bqual_ptr(bam_rec)[info.qpos] < 1) continue;
                            auto seq = bam_seq_ptr(bam_rec);
                            uint8_t base = bam_unpack_base(seq, info.qpos);
                            if (base == curr_vcf->ref) curr_vcf->ref_count[i_file]++;
                            else if (base == curr_vcf->alt) curr_vcf->alt_count[i_file]++;
                            else curr_vcf->err_count[i_file]++;
                        }
                    }
                }

                curr_vcf++;
            }

            return true;
        });
    });
    
    timestamp(clk_start, "pileup finished");

    std::stringstream outstream;

    for (auto& rec : records) {

        bool fails_min = !(rec.depth[0] >=0 && rec.depth[1] >= 0);
        bool is_not_zero = rec.ref_count[0] + rec.alt_count[0] + rec.err_count[0] +
                rec.ref_count[1] + rec.alt_count[1] + rec.err_count[1] > 0;
        
        if (is_not_zero && !fails_min) {
            if (bam_has_prefix) outstream << "chr";
            outstream << rec.chrom;
            outstream << "," << rec.pos + 1; // convert to VCF 1-base
            outstream << "," << byte2base_lo(rec.ref);
            outstream << "," << byte2base_lo(rec.alt);
            outstream << "," << rec.ref_count[0];
            outstream << "," << rec.alt_count[0];
            outstream << "," << rec.err_count[0];
            outstream << "," << rec.del_count[0];
            outstream << "," << rec.ref_count[1];
            outstream << "," << rec.alt_count[1];
            outstream << "," << rec.err_count[1];
            outstream << "," << rec.del_count[1];

            outstream<<std::endl;
        }
    }

    std::cout << outstream.rdbuf();

    timestamp(clk_start, "output written");

    index_free(index1);
    index_free(index2);

    return 0;

}

void load_vcf(std::vector<vcf_record>& records, std::vector<std::pair<size_t, size_t>>& work_ranges, const char *filename, const indices_t& indices, const std::vector<std::map<std::string, uint32_t>>& tid_maps) {
    
    uint8_t base2byte[20];
    base2byte['A'-'A'] = 1;
    base2byte['C'-'A'] = 2;
    base2byte['G'-'A'] = 4;
    base2byte['T'-'A'] = 8;

    // figure out vcf file size, allocate buffer, and read file
    struct stat stat_buf;
    int rc = stat(filename, &stat_buf);
    if(rc != 0) throw std::runtime_error("cannot get vcf file size");
    auto vcf_size = stat_buf.st_size;

    char *vcf_content = (char *)malloc(vcf_size);
    auto vcf_mf = mfile_open(filename);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, vcf_size), [&](tbb::blocked_range<size_t>& r) {
        memcpy(vcf_content + r.begin(), begin(vcf_mf) + r.begin(), r.end() - r.begin());
    });


    // sample 1% of the vcf buffer to estimate number of lines
    char *l_begin = vcf_content;
    char *l_end = strchr(vcf_content, '\n');

    size_t lines = 0;
    size_t sample_size = vcf_size / 100;
    for(; l_end != NULL && l_end - vcf_content < sample_size; l_begin = l_end + 1, l_end = strchr(l_begin, '\n')) lines++;
    lines *= 100;

    records.reserve(lines);

    // parse lines into vcf records
    l_begin = vcf_content;
    l_end = strchr(vcf_content, '\n');

    // -- START of data structure for breaking variants into work ranges -- //
    bool first = true;
    size_t batch_start = 0;
    size_t batch_end = 0;

    const size_t bam_span_limit = 1 * 1024 * 1024;
    work_ranges.reserve(lines / 330);

    const index_t& index0 = indices[0];
    const index_t& index1 = indices[1];

    uint32_t ref_id0 = 0;
    uint32_t ref_id1 = 0;

    uint64_t coffset_begin0 = 0;
    uint64_t coffset_end0= 0;

    uint64_t coffset_begin1 = 0;
    uint64_t coffset_end1= 0;

    // -- END of data structure for breaking variants into work ranges -- //

    char *chr, *pos, *id, *ref, *alt;
    for(;l_end != NULL; l_begin = l_end + 1, l_end = strchr(l_begin, '\n')) {

        if(*l_begin == '#') continue;
        chr = strchr(l_begin, '\t');
        pos = strchr(chr+1, '\t');
        id  = strchr(pos+1, '\t');
        ref = strchr(id+1, '\t');
        alt = strchr(ref+1, '\t');

        if(chr == NULL || pos == NULL || ref == NULL || alt == NULL) continue;
        if(ref - id > 2 || alt - ref > 2) continue;
        records.push_back({
                .chrom = std::string(l_begin, chr),
                .pos = (int32_t)strtol(chr+1, NULL, 0) - 1,
                .ref = base2byte[*(ref - 1) - 'A'],
                .alt = base2byte[*(alt - 1) - 'A']});

        // -- START of breaking variants into work ranges -- //
        if(first) {
            ref_id0 = tid_maps[0].at(records[0].chrom);
            ref_id1 = tid_maps[1].at(records[0].chrom);
            coffset_begin0 = index_coffset(index0.ref[ref_id0].ioffset[records[0].pos / 16384]);
            coffset_begin1 = index_coffset(index1.ref[ref_id1].ioffset[records[0].pos / 16384]);
            first = false;
        }
        else {
            coffset_end0 = index_coffset(index0.ref[ref_id0].ioffset[records[batch_end].pos / 16384]);
            coffset_end1 = index_coffset(index1.ref[ref_id1].ioffset[records[batch_end].pos / 16384]);

            if(records[batch_end].chrom != records[batch_start].chrom) {
                work_ranges.push_back({batch_start, batch_end});
                batch_start = batch_end;
                ref_id0 = tid_maps[0].at(records[batch_start].chrom);
                ref_id1 = tid_maps[1].at(records[batch_start].chrom);
                coffset_begin0 = coffset_end0;
                coffset_begin1 = coffset_end1;
            }
            else if(coffset_end0 - coffset_begin0 > bam_span_limit ||
                    coffset_end1 - coffset_begin1 > bam_span_limit) {
                work_ranges.push_back({batch_start, batch_end});
                batch_start = batch_end;
                ref_id0 = tid_maps[0].at(records[batch_start].chrom);
                ref_id1 = tid_maps[1].at(records[batch_start].chrom);
                coffset_begin0 = coffset_end0;
                coffset_begin1 = coffset_end1;
            }
        }
        // -- START of breaking variants into work ranges -- //

        batch_end++;
    }

    if (batch_start != batch_end)
        work_ranges.push_back({batch_start, batch_end});

    free(vcf_content);
}
