set(SUBGRAPH_MATCHING_FILE_LIST
        CommandList.cpp
        CommandList.h
        BuildTable.cpp
        BuildTable.h
        GraphQL.cpp
        GraphQL.h
        main.cpp)

add_executable(SubgraphMatching.out ${SUBGRAPH_MATCHING_FILE_LIST})
target_link_libraries(SubgraphMatching.out graph utility)
