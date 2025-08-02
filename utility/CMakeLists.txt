set(UTILITY_SRC
        Operation.cpp
        Operation.h
        CommandParser.cpp
        CommandParser.h
        ComputeSetIntersection.h
        ComputeSetIntersection.cpp
        analyze_symmetry.h
        analyze_symmetry.cpp
)

add_library(utility SHARED
        ${UTILITY_SRC})