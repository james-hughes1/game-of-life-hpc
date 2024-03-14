add_library(conway_lib
    STATIC
    world.cpp
    matrix.cpp
    parallel.cpp
)

target_include_directories(conway_lib
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

find_package(MPI REQUIRED)
target_link_libraries(conway_lib PUBLIC MPI::MPI_CXX)
