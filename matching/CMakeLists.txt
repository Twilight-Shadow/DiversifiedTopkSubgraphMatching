set(SUBGRAPH_MATCHING_FILE_LIST
        CommandList.cpp
        CommandList.h
        BuildTable.cpp
        BuildTable.h
        IntraMatch.cpp
        IntraMatch.h
        InterMatch.cpp
        InterMatch.h
        main.cpp)

add_executable(SubgraphMatching.out ${SUBGRAPH_MATCHING_FILE_LIST})
target_link_libraries(SubgraphMatching.out graph utility)
