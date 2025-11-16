set(HAN_SRC
    utils/util.cpp
    intersection_algos.cpp)

add_library(han SHARED ${HAN_SRC})

set_target_properties(han PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
