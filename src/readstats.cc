/***
 * (auto) parallel reads (coverage) count for some special regions
 * other reference genomic tools:
 * covtobed, samtools bedcov, bedtools multicov ...
 * Created by Schaudge King on 2021/10/18.
 */
#include <string>
#include <fstream>
#include <cstring>
#include <mmbam/mfile.h>
#include <mmbam/bam.h>
#include <mmbam/index.h>


#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
int main(int argc, const char *const argv[]) {

    if (argc < 2) {
        std::cout<<"usage: "<< argv[0] <<" <bam-file> [bai-file]"<<std::endl;
        exit(0);
    }

    auto bam_filename = argv[1];
    auto bai_filename = argc < 3 ? std::string(argv[1]).substr(0, strlen(argv[1]) - 4) + ".bai" : std::string(argv[2]);

    auto mfile = mfile_open(bam_filename);      // open BAM file as a mfile
    auto index = index_read(std::ifstream(bai_filename));  // open and parse index file

    size_t read_counts = 0;
    // create iteration intervals
    auto regions = index_to_regions(index, mfile->size);

    #pragma omp parallel for reduction(+:read_counts)
    for (size_t i=0; i < regions.size(); i++) {
        if(regions[i].first == regions[i].second) throw std::runtime_error("same region start and end");
        auto bam_records = bam_load_block(mfile, regions[i].first, regions[i].second);
        auto record_count = bam_count_records(bam_records);
        read_counts += record_count;
    }

    /***
     * uint32_t region_start = 0;
     * uint32_t region_end   = 55242466;

     * auto buffer = bam_load_region(mfile, index, 6, region_start, region_end);

     * // buffer now contains all reads on the loaded region.
     * bam_iterator bam_it(buffer);
     * bam_iterator bam_end(buffer, buffer.size());

     * while (bam_it < bam_end) {
     *   // process the read, get read name for example
     *   // std::cout << bam_read_name(bam_it) << std::endl;
     *   ++ read_counts;
     *   // advance the iterator
     *   bam_it++;
     *  }
     */

    std::cout << "Total reads counts in region is: " << read_counts << std::endl;

    index_free(index);
    return 0;
}
#pragma clang diagnostic pop
