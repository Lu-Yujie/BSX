add_subdirectory(a1)
add_subdirectory(a2)
add_subdirectory(a3)

add_library(aaa INTERFACE)

target_link_libraries(aaa INTERFACE a1 a2 a3)

set_target_properties(aaa PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
