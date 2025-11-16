set(ND_SRC
    bucket.cpp
    hierarchy.cpp
    k3-4.cpp
    kcore.cpp
    ktruss.cpp
    util.cpp
    nd_interface.cpp)

add_library(nd SHARED ${ND_SRC})
target_link_libraries(nd graph)

set_target_properties(nd PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
