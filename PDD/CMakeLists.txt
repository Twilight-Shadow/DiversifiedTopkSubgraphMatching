set (CMAKE_C_COMPILER "/home/hnu/Disk0/monica/hyc/gcc940/bin/gcc")
set (CMAKE_CXX_COMPILER "/home/hnu/Disk0/monica/hyc/gcc940/bin/g++")
# set (CMAKE_CXX_COMPILER "g++")
project(ParallelSubgraphMatching)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -g -Wall -march=native -mavx2 -pthread -mcmodel=large")

add_subdirectory(configuration)
add_subdirectory(graph)
add_subdirectory(matching)
add_subdirectory(utility)
add_subdirectory(ctpl)