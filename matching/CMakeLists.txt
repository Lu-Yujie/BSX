set(SUBGRAPH_MATCHING_FILE_LIST
        matchingcommand.cpp
        FilterVertices.cpp
        BuildEdgeIndex.cpp
        GenerateQueryPlan.cpp
        EvaluateQuery.cpp
        GenerateFilteringPlan.cpp
        StudyPerformance.cpp)

add_executable(BS ${SUBGRAPH_MATCHING_FILE_LIST})
target_link_libraries(BS graph utility)
target_link_libraries(BS gmp)
