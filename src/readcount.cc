#include<iostream>
#include<fstream>

#include "mmbam/mfile.h"
#include "mmbam/bam.h"
#include "mmbam/index.h"

// perform parallel read of a bam file.

#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"
int main(int argc, const char *const argv[]) {

    if(argc < 2) {
        std::cout<<"usage: "<<argv[0]<<" <bam-file> [bai-file]"<<std::endl;
        exit(0);
    }

    auto bam_filename = argv[1];
    auto bai_filename = argc < 3 ? std::string(argv[1]) + ".bai" : std::string(argv[2]);

    auto mfile = mfile_open(bam_filename);
    auto index = index_read(std::ifstream(bai_filename));

    uint64_t records = 0;

    // create iteration intervals
    auto regions = index_to_regions(index, mfile->size);

    #pragma omp parallel for reduction(+:records)
    for(size_t i=0; i<regions.size(); i++) {
        if(regions[i].first == regions[i].second) throw std::runtime_error("same region start and end");
        auto bam_records = bam_load_block(mfile, regions[i].first, regions[i].second);
        auto record_count = bam_count_records(bam_records);
        records += record_count;
    }

    std::cout<<"total records="<<records<<std::endl;

    index_free(index);
    return 0;
}
#pragma clang diagnostic pop
