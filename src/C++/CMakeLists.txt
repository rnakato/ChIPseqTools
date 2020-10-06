#cmake_minimum_required(VERSION 2.8)
project(ChIPseqTools CXX)

if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set( CMAKE_CXX_STANDARD 11 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )
set( CMAKE_CXX_EXTENSIONS OFF )
set( CMAKE_CXX_FLAGS "-std=c++11 -O2 -W -Wall")

include_directories("/usr/local/include")
link_directories("/usr/local/lib")

set(Boost_USE_STATIC_LIBS ON)
#set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

if(APPLE)
  set(BOOST_LIBS "-lpthread -lboost_thread-mt -lboost_program_options -lboost_system -lboost_filesystem")
else()
  set(BOOST_LIBS "-lpthread -lboost_thread -lboost_program_options -lboost_system -lboost_filesystem")
endif()

add_executable(compare_bed2loop compare_bed2loop.cpp gene_bed.cpp)
target_link_libraries(compare_bed2loop ssp_common common -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)

add_executable(gtf2refFlat gtf2refFlat.cpp)
target_link_libraries(gtf2refFlat common dd_func ssp_common -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)

add_executable(parseChIADropReadList parseChIADropReadList.cpp)
target_link_libraries(parseChIADropReadList common ssp_common -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)

add_executable(compare_bed2tss compare_bed2tss.cpp gene_bed.cpp)
target_link_libraries(compare_bed2tss common ssp_common dd_func -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)

add_executable(peak_occurance peak_occurance.cpp gene_bed.cpp statistics.cpp)
#target_link_libraries(peak_occurance common ssp_common dd_func -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)
target_link_libraries(peak_occurance common ssp_common dd_func ${BOOST_LIBS} -lgsl -lgslcblas)

#add_executable(multibed2gene multibed2gene.cpp gene_bed.cpp)
#target_link_libraries(compare_bed2tss common ssp_common dd_func -lboost_system -lboost_filesystem -lboost_iostreams -lboost_program_options -lgsl -lgslcblas)

#add_executable(FRiR FRiR.cpp)
#target_include_directories(FRiR
#  PUBLIC ${PROJECT_SOURCE_DIR}/DROMPAplus/submodules/SSP/common/
#)
#target_link_libraries(compare_bed2tss common ssp_common dd_func ssp_func -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)
