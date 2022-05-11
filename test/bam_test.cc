#include <gtest/gtest.h>
#include "mmbam/mfile.h"
#include "mmbam/mbgzf.h"
#include "mmbam/bam.h"
#include "mmbam/index.h"
#include <algorithm>
#include <fstream>
#include <stdio.h>

template<class T>
class error;

TEST(BAM, CanDetectCompleteHeader) {
    auto mfile = mfile_open("data/header_only.bam");
    bgzf_mfile_proxy_t bgzf_proxy(mfile);

    std::vector<uint8_t> bam_buffer;

    EXPECT_FALSE(bam_buffer_contains_header(bam_buffer));

    for(auto& bgzf_block : bgzf_proxy) {
        auto de_buff = bgzf_inflate(bgzf_block);
        bam_buffer.insert(bam_buffer.end(), de_buff.begin(), de_buff.end());
    }

    EXPECT_TRUE(bam_buffer_contains_header(bam_buffer));

    bam_buffer.resize(bam_buffer.size() / 2);
    bam_buffer.shrink_to_fit();

    EXPECT_FALSE(bam_buffer_contains_header(bam_buffer));
}

TEST(BAM, CanDecompressUsingPredicates) {
    auto mfile = mfile_open("data/header_only.bam");
    bgzf_mfile_proxy_t bgzf_proxy(mfile);
    std::vector<uint8_t> bam_buffer;

    auto it = bgzf_proxy.begin();
    while(!bam_buffer_contains_header(bam_buffer) && it != bgzf_proxy.end()) {
        auto de_buff = bgzf_inflate(*it);
        bam_buffer.insert(bam_buffer.end(), de_buff.begin(), de_buff.end());
        ++it;
    }

    EXPECT_TRUE(bam_buffer_contains_header(bam_buffer));
}

TEST(BAM, CanDecompressUsingTransform) {
    auto mfile = mfile_open("data/10blks.bam");
    bgzf_mfile_proxy_t bgzf_proxy(mfile);
    std::vector<uint8_t> bam_buffer;
    size_t header_blocks = 0;

    // decompress until end of header
    auto bam_rec_start = std::find_if(
            bgzf_proxy.begin(),
            bgzf_proxy.end(),
            [&bam_buffer, &header_blocks](auto& bgzf_block) {
                if(bam_buffer_contains_header(bam_buffer)) return true;
                auto de_buff = bgzf_inflate(bgzf_block);
                bam_buffer.insert(bam_buffer.end(), de_buff.begin(), de_buff.end());
                header_blocks++;
                return false;
            });

    auto first_bam_block = bgzf_inflate(*bam_rec_start);
    bam_rec_t* bam_rec = (bam_rec_t*)first_bam_block.data();

    EXPECT_EQ(bam_rec_start - bgzf_proxy.begin(), 26975);
    EXPECT_EQ(first_bam_block.size(), 65044);
    EXPECT_EQ(bam_rec->block_size, 188);
    EXPECT_EQ(bam_rec->ref_id, 0);
    EXPECT_EQ(bam_rec->pos, 10113);
}

TEST(BAM, CanDecompressOneInterval) {
    auto mfile = mfile_open("data/chr10.100blks.bam");
    auto index = index_read(std::ifstream("data/chr10.100blks.bam.bai"));

    auto bam_records = bam_load_block(mfile, index.ref[9].ioffset[61], index.ref[9].ioffset[62]);

    size_t total_reads = bam_count_records(bam_records);

    EXPECT_EQ(total_reads, 4535);
    index_free(index);
}

TEST(BAM, DecompressTwoIntervalEqualsTwoDecompression) {
    auto mfile = mfile_open("data/chr10.100blks.bam");
    auto index = index_read(std::ifstream("data/chr10.100blks.bam.bai"));


    auto bam_records1 = bam_load_block(mfile, index.ref[9].ioffset[61], index.ref[9].ioffset[63]);

    auto bam_records2_1 = bam_load_block(mfile, index.ref[9].ioffset[61], index.ref[9].ioffset[62]);
    auto bam_records2_2 = bam_load_block(mfile, index.ref[9].ioffset[62], index.ref[9].ioffset[63]);


    bam_records2_1.insert(bam_records2_1.end(), bam_records2_2.begin(), bam_records2_2.end());

    auto read_count1 = bam_count_records(bam_records1);
    auto read_count2 = bam_count_records(bam_records2_1);

    EXPECT_EQ(read_count1, 9299);
    EXPECT_EQ(read_count2, 9299);
}
    

TEST(BAMIterator, CanBeDereferenced) {
    auto mfile = mfile_open("data/10blks.bam");
    bgzf_mfile_proxy_t bgzf_proxy(mfile);
    std::vector<uint8_t> bam_buffer;

    auto bam_rec_start = std::find_if(
            bgzf_proxy.begin(),
            bgzf_proxy.end(),
            [&bam_buffer](auto& bgzf_block) {
                if(bam_buffer_contains_header(bam_buffer)) return true;
                auto de_buff = bgzf_inflate(bgzf_block);
                bam_buffer.insert(bam_buffer.end(), de_buff.begin(), de_buff.end());
                return false;
            });

    auto first_bam_block = bgzf_inflate(*bam_rec_start);
    bam_iterator bam_it(reinterpret_cast<bam_rec_t*>(first_bam_block.data()));

    EXPECT_EQ(bam_it->block_size, 188);
}

TEST(BAMIterator, CanBeAdvanced) {
    auto mfile = mfile_open("data/10blks.bam");
    bgzf_mfile_proxy_t bgzf_proxy(mfile);
    std::vector<uint8_t> bam_buffer;

    auto bam_rec_start = std::find_if(
            bgzf_proxy.begin(),
            bgzf_proxy.end(),
            [&bam_buffer](auto& bgzf_block) {
                if(bam_buffer_contains_header(bam_buffer)) return true;
                auto de_buff = bgzf_inflate(bgzf_block);
                bam_buffer.insert(bam_buffer.end(), de_buff.begin(), de_buff.end());
                return false;
            });

    auto first_bam_block = bgzf_inflate(*bam_rec_start);
    bam_iterator bam_it(reinterpret_cast<bam_rec_t*>(first_bam_block.data()));

    EXPECT_EQ(bam_it->block_size, 188);
    auto orig_it = bam_it;
    bam_it++;

    EXPECT_EQ(bam_it - orig_it, 192);
    EXPECT_EQ(bam_it->ref_id, 0);
    EXPECT_EQ(bam_it->pos, 10358);
}
