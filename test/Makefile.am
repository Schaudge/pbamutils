check_PROGRAMS = all_tests
TESTS = $(check_PROGRAMS)

AM_CPPFLAGS=-I$(top_srcdir)/src
AM_CXXFLAGS=$(GTEST_CFLAGS) $(PTHREAD_CFLAGS)
LDADD=$(GTEST_LIBS) ../src/libmmbam.la -ldeflate -ltbb $(PTHREAD_LIBS)

#mfile_test_SOURCES = mfile_test.cc
#mbgzf_test_SOURCES = mbgzf_test.cc
#bam_test_SOURCES = bam_test.cc
#index_test_SOURCES = index_test.cc

all_tests_SOURCES = mfile_test.cc mbgzf_test.cc index_test.cc bam_test.cc mpileup_test.cc
