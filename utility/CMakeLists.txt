add_subdirectory(bi_matching)
add_subdirectory(han)
add_subdirectory(nucleus_decomposition)

set(UTILITY_SRC
    graphoperations.cpp
    commandparser.cpp
    computesetintersection.cpp
    bitsetoperation.cpp)

add_library(utility SHARED ${UTILITY_SRC})
target_link_libraries(utility nd bigraphMatching han)

set_target_properties(utility PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
