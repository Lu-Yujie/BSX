set(VEQ_SRC
     veq.cpp)

add_library(VEQ SHARED ${VEQ_SRC})
target_link_libraries(VEQ graph)

set_target_properties(VEQ PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
)
