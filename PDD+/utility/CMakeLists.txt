set(UTILITY_SRC
        Operation.cpp
        Operation.h
        CommandParser.cpp
        CommandParser.h
        ComputeSetIntersection.h
        ComputeSetIntersection.cpp)

add_library(utility SHARED
        ${UTILITY_SRC})