cmake_minimum_required(VERSION 3.0)
project(SubgraphMatching)

include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/utility
    ${PROJECT_SOURCE_DIR}/sourceCode/VEQ
    ${PROJECT_SOURCE_DIR}/sourceCode/rapidMatch
)

set(CMAKE_CXX_FLAGS
    "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -g -Wall -march=native -pthread")
#     "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -g -Wall -march=native -fsanitize=address -pthread")

add_subdirectory(graph)
add_subdirectory(utility)
add_subdirectory(sourceCode)
add_subdirectory(matching)
