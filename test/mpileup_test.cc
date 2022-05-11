#include <gtest/gtest.h>
#include "mmbam/mfile.h"
#include "mmbam/mbgzf.h"
#include "mmbam/bam.h"
#include "mmbam/index.h"
#include "mmbam/mpileup.h"
#include <algorithm>
#include <fstream>
#include <stdio.h>

constexpr uint16_t flag_fail = READ_FAILED_QUALITY_CHECKS | READ_SECONDARY_ALIGNMENT | READ_PCR_OPTICAL_DUPLICATE | READ_UNMAPPED;
static bool true_predicate(const bam_rec_t&) { return true; }

TEST(MPILEUP_1FILE, CanPileupAtGivenLocation) {
    auto mfile = mfile_open("data/chr10.100blks.bam");
    auto index = index_read(std::ifstream("data/chr10.100blks.bam.bai"));

    mfiles_t mfiles{mfile};
    indices_t indices{index};

    size_t coverage = 0;

    // chr10, pos: 1001280
    mpileup(mfiles, indices, 9, 1001279, 1001280, true_predicate, [&coverage](auto &p) {
            if(p.pos > 1001279) return false;
            if(p.pos == 1001279) coverage = p.depth[0];
            return true;
            });

    EXPECT_EQ(coverage, 88); // as calculated from samtools view

}

TEST(MPILEUP_1FILE, CanPileupAtEmptyLocation) {
    auto mfile = mfile_open("data/chr10.100blks.bam");
    auto index = index_read(std::ifstream("data/chr10.100blks.bam.bai"));
    
    mfiles_t mfiles{mfile};
    indices_t indices{index};

    size_t coverage = 0;

    // chr10, pos: 14101
    mpileup(mfiles, indices, 9, 14100, 14101, true_predicate, [&coverage](auto &p) {
            if(p.pos > 14100) return false;
            if(p.pos == 14100) coverage = p.depth[0];
            return true;
            });

    EXPECT_EQ(coverage, 0);
}

TEST(MPILEUP_1FILE, CanCalculateQueryPosition) {
    auto mfile = mfile_open("data/chr10.100blks.bam");
    auto index = index_read(std::ifstream("data/chr10.100blks.bam.bai"));

    mfiles_t mfiles{mfile};
    indices_t indices{index};

    // chr10, pos: 1001280
    mpileup(mfiles, indices, 9, 1001279, 1001280, true_predicate, [](auto &p) {
            if(p.pos > 1001279) return false;
            if(p.pos == 1001279) {
                bam_iterator bam_it(*(p.reads_buffer[0]));
                bam_iterator bam_end(*(p.reads_buffer[0]), p.reads_buffer[0]->size());

                // first read should have a query position of 123
                EXPECT_EQ(p.get_info(0, 0).qpos, 248);

                for(size_t i=0; i<p.depth[0]; i++) {
                    auto seq = bam_seq_ptr(BAMREF(p.reads_buffer[0]->data() + p.get_info(0, i).offset));
                    uint8_t base = bam_unpack_base(seq, p.get_info(0, i).qpos);
                    EXPECT_EQ(base, 4); //G
                }
            }
            return true;
            });
}

TEST(MPILEUP_1FILE, CanPileupHeterozygousPosition) {
    auto mfile = mfile_open("data/chr10.100blks.bam");
    auto index = index_read(std::ifstream("data/chr10.100blks.bam.bai"));
    
    mfiles_t mfiles{mfile};
    indices_t indices{index};

    //chr10, pos: 15065
    mpileup(mfiles, indices, 9, 15064, 15065, true_predicate, [](auto &p) {
        if(p.pos > 15064) return false;
        if(p.pos == 15064) {
            // coverage without filtering
            EXPECT_EQ(p.depth[0], 64);

            int base_counter[16];
            memset(base_counter, 0, sizeof(int) * 16);

            for(size_t read_idx = 0; read_idx < p.depth[0]; read_idx++) {
                const bam_rec_t* bam_rec = reinterpret_cast<const bam_rec_t *>(
                        p.reads_buffer[0]->data() + p.get_info(0, read_idx).offset);

                if(bam_rec->flag & flag_fail) {
                    continue;
                }

                auto seq = bam_seq_ptr(bam_rec);
                uint8_t base = bam_unpack_base(seq, p.get_info(0, read_idx).qpos);

                base_counter[base]++;
            }

            // with samtools default filter, 14 As and 15 Gs
            EXPECT_EQ(base_counter[1], 21); // 14 As
            EXPECT_EQ(base_counter[4], 43); // 15 Gs
        }
        return true;
    });
}

TEST(MPILEUP_1FILE, CanPileupInDeletionRegions) {
    auto mfile = mfile_open("data/chr10.100blks.bam");
    auto index = index_read(std::ifstream("data/chr10.100blks.bam.bai"));

    mfiles_t mfiles{mfile};
    indices_t indices{index};

    // chr10, pos: 417438; 31 reads contain deleted reference seq
    mpileup(mfiles, indices, 9, 417437, 417438, true_predicate, [](auto &p) {
        if(p.pos > 417437) return false;
        if(p.pos == 417437) {
            int is_del_count = 0;

            // coverage without filtering
            EXPECT_EQ(p.depth[0], 51);

            for(size_t i=0; i<p.depth[0]; i++)
                is_del_count += p.get_info(0, i).is_deletion ? 1 : 0;

            EXPECT_EQ(is_del_count, 31);
        }
        return true;
    });
}

TEST(MPILEUP_2FILES, CanPileupHeterozygousPosition) {
    auto mfile1 = mfile_open("data/chr10.100blks.bam");
    auto mfile2 = mfile_open("data/chr10.100blks.2.bam");

    auto index1 = index_read(std::ifstream("data/chr10.100blks.bam.bai"));
    auto index2 = index_read(std::ifstream("data/chr10.100blks.2.bam.bai"));

    mfiles_t mfiles{mfile1, mfile2};
    indices_t indices{index1, index2};

    //chr10, pos: 15065
    mpileup(mfiles, indices, 9, 15064, 15065, true_predicate, [](auto &p) {
        if(p.pos > 15064) return false;
        if(p.pos == 15064) {
            // coverage without filtering
            EXPECT_EQ(p.depth[0], 64);
            EXPECT_EQ(p.depth[1], 85);

            int base_counter[2][16];
            memset(base_counter, 0, sizeof(int) * 2 * 16);

            for(size_t file_idx = 0; file_idx < 2; file_idx++) {
                for(size_t read_idx = 0; read_idx < p.depth[file_idx]; read_idx++) {
                    const bam_rec_t* bam_rec = reinterpret_cast<const bam_rec_t *>(
                            p.reads_buffer[file_idx]->data() + p.get_info(file_idx, read_idx).offset);

                    if(bam_rec->flag & flag_fail) {
                        continue;
                    }

                    auto seq = bam_seq_ptr(bam_rec);
                    uint8_t base = bam_unpack_base(seq, p.get_info(file_idx, read_idx).qpos);

                    base_counter[file_idx][base]++;
                }
            }

            // with samtools default filter and -x, 14 As and 15 Gs
            EXPECT_EQ(base_counter[0][1], 21); // 14 As
            EXPECT_EQ(base_counter[0][4], 43); // 15 Gs
            EXPECT_EQ(base_counter[1][1], 44); // 14 As
            EXPECT_EQ(base_counter[1][4], 40); // 15 Gs

        }
        return true;
    });


}

TEST(MPILEUP_2FILES, CanPileupAtEmptyLocation) {
    auto mfile1 = mfile_open("data/chr10.100blks.bam");
    auto mfile2 = mfile_open("data/chr10.100blks.2.bam");

    auto index1 = index_read(std::ifstream("data/chr10.100blks.bam.bai"));
    auto index2 = index_read(std::ifstream("data/chr10.100blks.2.bam.bai"));

    mfiles_t mfiles{mfile1, mfile2};
    indices_t indices{index1, index2};

    int coverage1 = 0;
    int coverage2 = 0;

    // chr10, pos: 500
    mpileup(mfiles, indices, 9, 499, 500, true_predicate, [&](auto &p) {
        if(p.pos > 499) return false;
        if(p.pos == 499) {
            // coverage without filtering
            coverage1 = p.depth[0];
            coverage2 = p.depth[1];
        }
        return true;
    });

    EXPECT_EQ(coverage1, 0);
    EXPECT_EQ(coverage2, 0);

}

TEST(MPILEUP_2FILES, CanPileupInDeletionRegions) {
    auto mfile1 = mfile_open("data/chr10.100blks.bam");
    auto mfile2 = mfile_open("data/chr10.100blks.2.bam");

    auto index1 = index_read(std::ifstream("data/chr10.100blks.bam.bai"));
    auto index2 = index_read(std::ifstream("data/chr10.100blks.2.bam.bai"));

    mfiles_t mfiles{mfile1, mfile2};
    indices_t indices{index1, index2};

    // chr10, pos: 417438
    mpileup(mfiles, indices, 9, 417437, 417438, true_predicate, [](auto &p) {
        if(p.pos > 417437) return false;
        if(p.pos == 417437) {
            int is_del_count = 0;

            // coverage without filtering
            EXPECT_EQ(p.depth[0], 51);
            EXPECT_EQ(p.depth[1], 55);

            for(size_t i=0; i<p.depth[0]; i++)
                is_del_count += p.get_info(0, i).is_deletion ? 1 : 0;

            EXPECT_EQ(is_del_count, 31); //without flag filter

            is_del_count = 0;
            for(size_t i=0; i<p.depth[1]; i++)
                is_del_count += p.get_info(1, i).is_deletion ? 1 : 0;

            EXPECT_EQ(is_del_count, 25); //without flag filter

        }
        return true;
    });
}
