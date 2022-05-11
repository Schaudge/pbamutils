#include <gtest/gtest.h>
#include <iostream>
#include <algorithm>

#include "mmbam/mfile.h"
#include "mmbam/mbgzf.h"

#include <tbb/parallel_for.h>

TEST(MBGZFIterator, CanCountBlocks) {
    auto mfile = mfile_open("data/10blks.bam");
    auto bgzf_proxy = bgzf_mfile_proxy_t(mfile);
    size_t blocks = 0;
    bool last_is_eof = false;

    auto it_begin = bgzf_proxy.begin();
    auto it_end   = bgzf_proxy.end();

    for(; it_begin != it_end; ++it_begin) {
        EXPECT_EQ(it_begin->id1, 31);
        EXPECT_EQ(it_begin->id2, 139);
        blocks++;
        last_is_eof = is_bgzf_eof_block(*it_begin);
    }

    EXPECT_EQ(blocks, 105);
    EXPECT_TRUE(last_is_eof);
}

TEST(MBGZFIterator, CanCountBlocksWithRange) {
    auto mfile = mfile_open("data/10blks.bam");
    auto bgzf_proxy = bgzf_mfile_proxy_t(mfile);
    size_t blocks = 0;
    bool last_is_eof = false;
    for(auto& bgzf_block : bgzf_proxy) {
        EXPECT_EQ(bgzf_block.id1, 31);
        EXPECT_EQ(bgzf_block.id2, 139);
        blocks++;
        last_is_eof = is_bgzf_eof_block(bgzf_block);
    }
    EXPECT_EQ(blocks, 105);
    EXPECT_TRUE(last_is_eof);
}

TEST(MBGZF, CanInflate) {
    auto mfile = mfile_open("data/10blks.bam");
    auto expected_size = 6701624;

    decltype(expected_size) total_size = 0;
    auto bgzf_proxy = bgzf_mfile_proxy_t(mfile);

    for(auto& bgzf_block : bgzf_proxy) {
        auto inflate_results = bgzf_inflate(bgzf_block);
        if(!is_bgzf_eof_block(bgzf_block)) EXPECT_GT(inflate_results.size(), 0);
        total_size += inflate_results.size();
    }
    EXPECT_EQ(expected_size, total_size);
}

TEST(MBGZF, CanParallelInflateSequentially) {
    auto mfile = mfile_open("data/10blks.bam");
    auto bgzf_proxy = bgzf_mfile_proxy_t(mfile);

    std::vector<uint8_t> ref_buffer;
    for(auto& bgzf_block : bgzf_proxy) {
        auto inflate_results = bgzf_inflate(bgzf_block);
        ref_buffer.insert(ref_buffer.end(), inflate_results.cbegin(), inflate_results.cend());
    }

    auto buffer = bgzf_inflate_range_p(bgzf_proxy.begin().ptr, mfile->size, [](
        auto* src, auto src_len, auto dest_len,
        auto& src_off_vector, auto& dest_off_vector,
        auto inflate) -> std::vector<uint8_t> {

        std::vector<uint8_t> buffer;
        buffer.resize(dest_len);

        for(size_t i=0; i<src_off_vector.size(); i++) {
            bgzf_iterator_t block((bgzf_block_t*)(src + src_off_vector[i]));
            inflate(src + src_off_vector[i], block->bsize + 1,
                    buffer.data() + dest_off_vector[i], bgzf_isize(*block));
        }

        return buffer;
    });

    EXPECT_EQ(buffer.size(), ref_buffer.size());
    auto equal = std::equal(buffer.cbegin(), buffer.cend(), ref_buffer.cbegin());
    EXPECT_EQ(equal, true);
}

TEST(MBGZF, CanParallelInflate) {
    auto mfile = mfile_open("data/10blks.bam");
    auto bgzf_proxy = bgzf_mfile_proxy_t(mfile);

    std::vector<uint8_t> ref_buffer;
    for(auto& bgzf_block : bgzf_proxy) {
        auto inflate_results = bgzf_inflate(bgzf_block);
        ref_buffer.insert(ref_buffer.end(), inflate_results.cbegin(), inflate_results.cend());
    }

    auto buffer = bgzf_inflate_range_p(bgzf_proxy.begin().ptr, mfile->size, [](
        auto* src, auto src_len, auto dest_len,
        auto& src_off_vector, auto& dest_off_vector,
        auto inflate) -> std::vector<uint8_t> {

        std::vector<uint8_t> buffer;
        buffer.resize(dest_len);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, src_off_vector.size()), [&](const auto& r){
            for(size_t i=r.begin(); i<r.end(); i++) {
                bgzf_iterator_t block((bgzf_block_t*)(src + src_off_vector[i]));
                inflate(src + src_off_vector[i], block->bsize + 1,
                        buffer.data() + dest_off_vector[i], bgzf_isize(*block));
            }
        });

        return buffer;
    });

    EXPECT_EQ(buffer.size(), ref_buffer.size());
    auto equal = std::equal(buffer.cbegin(), buffer.cend(), ref_buffer.cbegin());
    EXPECT_EQ(equal, true);
}
