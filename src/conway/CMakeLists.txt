add_library(conway_lib
    STATIC
    conway.cpp
    matrix.cpp
)

target_include_directories(conway_lib
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)
