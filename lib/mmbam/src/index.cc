#include <stdint.h>
#include <fstream>
#include "mmbam/index.h"
#include "mmbam/bam.h"
#include "mmbam/mbgzf.h"
#include <vector>
#include <exception>

void index_free(index_t& index) {
    for(size_t i=0; i<index.n_ref; i++) {
        for(size_t j=0; j<index.ref[i].n_bin; j++) delete [] index.ref[i].bin[j].chunk;
        delete [] index.ref[i].bin;
        delete [] index.ref[i].ioffset;
    }
    delete [] index.ref;
}

std::vector<region> index_to_regions(const index_t& index, uint64_t filesize) {
    std::vector<region> regions;
    using ref_idx_t = decltype(index.n_ref);
    using intv_idx_t = decltype(index.ref[0].n_intv);

    bool first = true;
    uint64_t block_start = 0;
    uint64_t block_end = 0;

    for(ref_idx_t i_ref = 0; i_ref < index.n_ref; ++i_ref) {
        for(intv_idx_t i_intv = 0; i_intv < index.ref[i_ref].n_intv; ++i_intv) {
            if(first) {
                block_start = index.ref[i_ref].ioffset[i_intv];
                first = false;
            }

            block_end = index.ref[i_ref].ioffset[i_intv];
            if(block_end != block_start) {
                regions.push_back({block_start, block_end});
                block_start = block_end;
            }
        }
    }


    if(index_coffset(block_start) < filesize) {
        block_end = filesize << 16;
        regions.push_back({block_start, block_end});
    }

    return regions;
}
