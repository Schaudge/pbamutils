#include <gtest/gtest.h>
#include <fstream>
#include <iostream>
#include "mmbam/index.h"
#include "mmbam/bam.h"

TEST(INDEX, CanBeLoaded) {
    auto bai_stream = std::ifstream("data/chr10.100blks.bam.bai");
    auto index = index_read(bai_stream);
    EXPECT_EQ(index.magic[0], 'B');
    EXPECT_EQ(index.magic[1], 'A');
    EXPECT_EQ(index.magic[2], 'I');
    EXPECT_EQ(index.magic[3], '\1');
    EXPECT_GT(index.n_ref, 0);
    for(int i=0; i<index.n_ref; i++)
        if(i == 9) EXPECT_GT(index.ref[i].n_intv, 0);
        else       EXPECT_EQ(index.ref[i].n_intv, 0);
}

// test getting reads in a region
TEST(INDEX, GetsAllReadsInARegion) {

    // region to test: chr10:1,500,000-1,530,000
    // spans at least 2 linear blocks

    auto bai_stream = std::ifstream("data/chr10.100blks.bam.bai");
    auto index = index_read(bai_stream);

    auto mfile = mfile_open("data/chr10.100blks.bam");

    using intv_idx_t = decltype(index.ref[9].n_intv);

    uint32_t region_start = 1500000;
    uint32_t region_end   = 1530000;

    auto bam_buffer = bam_load_region(mfile, index, 9, region_start, region_end);
    bam_iterator bam_it(bam_buffer);
    bam_iterator bam_end(bam_buffer, bam_buffer.size());

    int hits = 0;

    // expect all reads to be within the specified region
    while(bam_it < bam_end) {
        EXPECT_TRUE(bam_it->pos < region_end);
        EXPECT_TRUE(bam_it->pos + bam_query_length(bam_it) >= region_start);
        hits++;
        bam_it++;
    }

    EXPECT_EQ(hits, 8676); // as calculated from samtools view
}

// test getting reads that overlaps a position
TEST(INDEX, GetsAllReadsOverAPosition) {
    // position to test: chr10:1001280
    //
    auto bai_stream = std::ifstream("data/chr10.100blks.bam.bai");
    auto index = index_read(bai_stream);
    auto mfile = mfile_open("data/chr10.100blks.bam");

    int32_t region_start = 1001280;
    int32_t region_end   = 1001281;

    auto bam_buffer = bam_load_region(mfile, index, 9, region_start, region_end);
    bam_iterator bam_it(bam_buffer);
    bam_iterator bam_end(bam_buffer, bam_buffer.size());

    // all reads should overlap chr10:1001280
    int hits = 0;
    while(bam_it < bam_end) {
        if(bam_it->flag & 0x4) {
            bam_it++;
            continue;
        }
        EXPECT_TRUE(bam_it->pos < region_end);
        EXPECT_TRUE(bam_it->pos + bam_query_length(bam_it) >= region_start);
        hits++;
        bam_it++;
    }
    EXPECT_EQ(hits, 87);
}
