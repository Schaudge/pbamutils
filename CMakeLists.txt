cmake_minimum_required(VERSION 2.8.8)
set(PROJECT_NAME_STR pbamutils)

project(${PROJECT_NAME_STR} C CXX)
add_compile_options(-Wall -std=c++14 -fopenmp)

# required header and source file
include_directories(${PROJECT_SOURCE_DIR}/lib/mmbam/src)

add_library(mmbam STATIC lib/mmbam/src/bam.cc lib/mmbam/src/index.cc 
	                 lib/mmbam/src/mbgzf.cc lib/mmbam/src/mpileup.cc)

add_executable(facets-snp-pileup src/snp-pileup.cc)
target_link_libraries(facets-snp-pileup -fopenmp mmbam -ldeflate -lz -ltbb)

add_executable(flagstats src/flagstats.cc)
target_link_libraries(flagstats -fopenmp mmbam -ldeflate -lz -ltbb)

add_executable(readcount src/readcount.cc)
target_link_libraries(readcount -fopenmp mmbam -ldeflate -lz -ltbb)

add_executable(readstats  src/readstats.cc)
target_link_libraries(readstats -fopenmp mmbam -ldeflate -lz -ltbb)

add_executable(pmpileup  src/mpileup_parallel.cc)
target_link_libraries(pmpileup -fopenmp mmbam -ldeflate -lz -ltbb)

#[test]
enable_testing()
find_package(GTest REQUIRED)
include(GoogleTest)

add_executable(tests test/index_test.cc test/bam_test.cc test/mfile_test.cc test/mpileup_test.cc)
target_link_libraries(tests GTest::GTest GTest::Main -fopenmp mmbam -ldeflate -lz -ltbb)
gtest_discover_tests(tests)

