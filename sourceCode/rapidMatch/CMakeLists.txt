# construct the lib of RapidMatch
add_subdirectory(nucleus_decomposition)

set(RM_SRC
    leapfrogtriejoin/leapfrogtriejoin.cpp
    hashjoin/nop_join.cpp
    primitive/search.cpp primitive/scan.cpp primitive/semi_join.cpp primitive/projection.cpp primitive/nlf_filter.cpp
    relation/catalog.cpp
    execution_tree/execution_tree_generator.cpp execution_tree/execution_tree.cpp execution_tree/execution_tree_node.cpp)

# add include path
include_directories(
    ${PROJECT_SOURCE_DIR}/sourceCode/rapidMatch
)

add_library(RM SHARED ${RM_SRC})

target_link_libraries(RM utility nd)

set_target_properties(RM PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
