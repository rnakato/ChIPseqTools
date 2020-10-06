#cmake_minimum_required(VERSION 2.8)
#project(ChIPseqTools CXX)

#set( CMAKE_CXX_STANDARD 14 )
#set( CMAKE_CXX_STANDARD_REQUIRED ON )
#set( CMAKE_CXX_EXTENSIONS OFF )
#set( CMAKE_CXX_FLAGS "-std=c++14 -O2 -W -Wall")

add_executable(compare_bed2loop compare_bed2loop.cpp gene_bed.cpp)
target_link_libraries(compare_bed2loop ssp_common common -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)

add_executable(gtf2refFlat gtf2refFlat.cpp)
target_link_libraries(gtf2refFlat common dd_func ssp_common -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)

add_executable(parseChIADropReadList parseChIADropReadList.cpp)
target_link_libraries(parseChIADropReadList common ssp_common -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)

add_executable(compare_bed2tss compare_bed2tss.cpp gene_bed.cpp)
target_link_libraries(compare_bed2tss common ssp_common dd_func -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)

add_executable(peak_occurance peak_occurance.cpp gene_bed.cpp statistics.cpp)
target_link_libraries(peak_occurance common ssp_common dd_func -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)

#add_executable(multibed2gene multibed2gene.cpp gene_bed.cpp)
#target_link_libraries(compare_bed2tss common ssp_common dd_func -lboost_system -lboost_filesystem -lboost_iostreams -lboost_program_options -lgsl -lgslcblas)

#add_executable(FRiR FRiR.cpp)
#target_include_directories(FRiR
#  PUBLIC ${PROJECT_SOURCE_DIR}/DROMPAplus/submodules/SSP/common/
#)
#target_link_libraries(compare_bed2tss common ssp_common dd_func ssp_func -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams -lgsl -lgslcblas)
