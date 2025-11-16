cmake_minimum_required(VERSION 3.0)
project(MyProject)

include_directories(
    ${CMAKE_SOURCE_DIR}
)

# add subdirectory a & b
add_subdirectory(a)
add_subdirectory(b)

# source file
set(SOURCES main.cpp)

add_executable(main ${SOURCES})

# link all sub-lib to main
target_link_libraries(main aaa bbbb)

# output path
set_target_properties(main PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)
