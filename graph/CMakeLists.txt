set(GRAPH_SRC
    graph.cpp)

add_library(graph SHARED
    ${GRAPH_SRC})

# output path
set_target_properties(graph PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
