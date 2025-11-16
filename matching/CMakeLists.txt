set(SUBGRAPH_MATCHING_FILE_LIST
    matchingcommand.cpp
    FilterVertices.cpp
    BuildEdgeIndex.cpp
    GenerateQueryPlan.cpp
    EvaluateQuery.cpp
    GenerateFilteringPlan.cpp
    Experiment.cpp)

include_directories(
    ${PROJECT_SOURCE_DIR}/sourceCode
)

add_executable(BS ${SUBGRAPH_MATCHING_FILE_LIST})
target_link_libraries(BS RM VEQ)
target_link_libraries(BS gmp)

# set ouput path
set_target_properties(BS PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
