cmake_minimum_required(VERSION 2.8)
project(ChIPseqTools CXX)

#set( CMAKE_CXX_STANDARD 14 )
#set( CMAKE_CXX_STANDARD_REQUIRED ON )
#set( CMAKE_CXX_EXTENSIONS OFF )
#set( CMAKE_CXX_FLAGS "-std=c++14 -O2 -W -Wall")

include_directories("/usr/local/include")

add_subdirectory(src/C)
add_subdirectory(submodules/DROMPAplus)
add_subdirectory(src/C++)
