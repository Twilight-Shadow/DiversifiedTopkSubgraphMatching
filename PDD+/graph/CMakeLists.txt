set(GRAPH_SRC
        BasicGraph.cpp
        MatchGraph.cpp)

add_library(graph SHARED
        ${GRAPH_SRC})