/***
 * (fast) basic mpileup by position (parallel) use mmbam library
 * Written By Schaudge King
 * 2022-03-09
 */
#include <fstream>
#include <tbb/parallel_for_each.h>
#include <mmbam/mfile.h>
#include <mmbam/bam.h>
#include <mmbam/index.h>
#include <mmbam/mpileup.h>

struct work_items {
    // location to pileup
    int32_t ref_id, pos;
    // the reference base at this location
    uint8_t ref_base;

    // results after processing the piled up location
    // This is necessary because parallel jobs will need independent
    // locations to write their results, else a mutex will be necessary
    int ref_count[2] = {0, 0};     // count ref alleles per input buffer
    int other_count[2] = {0, 0 };  // count other alleles per input buffer
};

int main(int argc, const char *const argv[]) {
    // mpileup engine supports multiple input files.
    // This example shows how to pileup two files simultaneously

    auto mfile1 = mfile_open("/home/yuanshenran/datasets/bam/gzy201001024_g_13_tumor_recal.bam");
    auto mfile2 = mfile_open("/home/yuanshenran/datasets/bam/gzy201002011_g_16_tumor_recal.bam");

    auto index1 = index_read(std::ifstream("/home/yuanshenran/datasets/bam/gzy201001024_g_13_tumor_recal.bai"));
    auto index2 = index_read(std::ifstream("/home/yuanshenran/datasets/bam/gzy201002011_g_16_tumor_recal.bai"));

    mfiles_t mfiles{mfile1, mfile2};
    indices_t indices{index1, index2};

    std::vector<work_items> parallel_items = {
            {ref_id: 2, pos: 178936091, ref_base: 0x2},
            {ref_id: 6, pos: 55249071, ref_base: 0x1}
    };
    // the following line generate the (const) parallel loop index,
    // which used to fill the ref_id, pos, and ref_base fields of each work_items,
    // and initialize the ref_count and other_count fields to 0
    std::vector<int> parallel_index{0,1};

    tbb::parallel_for_each(
        parallel_index.cbegin(),
        parallel_index.cend(),
        [&](auto& idx) {

            // the filter lambda returns false if the mapping quality of
            // a read is below a given a hardcoded threshold (1)
            auto filter_predicate = [](const auto& bam_rec) {
                if (bam_rec.mapq < 1) return false;
                return true;
            };

            auto r = parallel_items.begin() + idx;
            // The visitor function lambda will be called at each piled up
            // location with a mpileup_t struct as the parameter
            auto visitor_func = [&r](const auto& p) {

                // skip all positions before the desired position
                if (p.pos < r->pos) return true;

                // halt pileup if we are past the desired position
                if (p.pos > r->pos) return false;

                // here implies that the piled up position is equal to r.pos
                // count reference and other reads per input buffer
                for (size_t i_file = 0; i_file < 2; i_file++) {
                    auto* buffer = p.reads_buffer[i_file]->data();

                    // iterate over reads in buffer i_file
                    for (size_t i = 0; i < p.depth[i_file]; i++) {
                        auto& info = p.get_info(i_file, i);
                        const bam_rec_t* bam_record = BAMREF(buffer + info.offset);

                        // only count non-deletion alleles
                        if (!info.is_deletion) {
                            if (bam_bqual_ptr(bam_record)[info.qpos] < 1)
                                continue; // ignore low quality bases
                            auto seq = bam_seq_ptr(bam_record);
                            auto base = bam_unpack_base(seq, info.qpos);
                            if (base == r->ref_base) r->ref_count[i_file]++;
                            else r->other_count[i_file]++;
                        }
                    } // end for-each read
                } // end for-each buffer

                return true;
            };

            // call mpileup engine
            mpileup(mfiles, indices, r->ref_id, r->pos, r->pos+1,
                    filter_predicate, visitor_func);

    });

    // parallel_items now contain the results
    // post-processing
    for(const auto& item: parallel_items) {
        std::cout << "ref: " << item.ref_count[0] << ", alt: " << item.other_count[0] << "\n";
    }

    index_free(index1);
    index_free(index2);

    return 0;
}

