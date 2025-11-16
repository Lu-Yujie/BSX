add_library(bigraphMatching SHARED
    bigraphMatching.cpp
)

set_target_properties(bigraphMatching PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
