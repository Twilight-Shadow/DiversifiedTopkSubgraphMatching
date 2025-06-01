project(ParallelSubgraphMatching)

include_directories(
        ${PROJECT_SOURCE_DIR}/
)

set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++17 -O3 -g -Wall -march=native -pthread -mcmodel=large")

add_subdirectory(configuration)
add_subdirectory(graph)
add_subdirectory(matching)
add_subdirectory(utility)
add_subdirectory(ctpl)