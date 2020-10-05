cmake_minimum_required(VERSION 2.8)
project(ChIPseqTools C)

set( CMAKE_C_FLAGS "-Wall -g -W -O3 -lm")

include_directories("/usr/local/include")

add_executable(compare_bs compare_bs.c compare.c readgene.c my.c stringp.c)
target_link_libraries(compare_bs m)

add_executable(compare_bs2CRM compare_bs.c compare.c readgene.c my.c stringp.c)
target_link_libraries(compare_bs2CRM m)
