#include<gtest/gtest.h>
#include<string.h>
#include "mmbam/mfile.h"

constexpr auto test_bam = "data/10blks.bam";

TEST(MFile, CanOpenFile) {
    auto mfile = mfile_open(test_bam);
    EXPECT_NE(mfile.get(), nullptr);
    EXPECT_NE(mfile->mmptr, nullptr);
    EXPECT_NE(mfile->fd, -1);
    EXPECT_EQ(mfile->size, 2082899);
}

TEST(MFile, AutoCloseFile) {
    int fd;
    {
        auto mfile = mfile_open(test_bam);
        fd = mfile->fd;
    }
    auto checkval = fcntl(fd, F_GETFD);
    EXPECT_TRUE(checkval == -1 || errno == EBADF);
}

TEST(MFileIterator, CanReturn) {
    auto mfile = mfile_open(test_bam);
    auto *m_begin = begin<uint8_t>(mfile);
    EXPECT_EQ(*m_begin++, 31);      // ID1
    EXPECT_EQ(*m_begin++, 139);     // ID2
    EXPECT_EQ(*m_begin++, 8);       // CM
    EXPECT_EQ(*m_begin++, 4);       // FLG
}

TEST(MFileIterator, CanIterateOverEntireFile) {
    auto mfile = mfile_open(test_bam);
    off_t size = 0;
    for(auto& byte : mfile) {
        size++;
    }
    EXPECT_EQ(size, mfile->size);
}

TEST(MFileIterator, CanDetectBGZFEof) {
    auto mfile = mfile_open(test_bam);
    auto *m_eofmk = end<uint8_t>(mfile);
    uint8_t eof_bytes[] = {
        0x1f, 0x8b, 0x08, 0x04,
        0x00, 0x00, 0x00, 0x00,
        0x00, 0xff, 0x06, 0x00,
        0x42, 0x43, 0x02, 0x00,
        0x1b, 0x00, 0x03, 0x00,
        0,0,0,0,0,0,0,0
    };
    m_eofmk -= 28;
    EXPECT_EQ(strncmp((const char*)m_eofmk, (const char*)eof_bytes, 28), 0);
}
