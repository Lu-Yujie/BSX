add_library(bbbb SHARED
    b1.cpp
    b2.cpp
)

target_link_libraries(bbbb aaa)

set_target_properties(bbbb PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
