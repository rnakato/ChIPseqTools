cmake_minimum_required(VERSION 2.8)
project(ChIPseqTools CXX)

if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
      cmake_policy(SET CMP0015 NEW)
endif(COMMAND cmake_policy)

add_subdirectory(src/C)
add_subdirectory(submodules/DROMPAplus)
add_subdirectory(src/C++)
