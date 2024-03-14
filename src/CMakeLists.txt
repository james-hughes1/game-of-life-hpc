add_subdirectory(conway)
add_subdirectory(timing)

add_executable(main main.cpp)
ClangFormat(main .)

add_executable(time_simulation time_simulation.cpp)
ClangFormat(time_simulation .)

add_executable(time_count_neighbours time_count_neighbours.cpp)
ClangFormat(time_count_neighbours .)

add_executable(parallel parallel.cpp)
ClangFormat(parallel .)

add_executable(time_1d_decomp
time_1d_decomp.cpp
)
ClangFormat(time_1d_decomp .)

target_link_libraries(main PRIVATE conway_lib)
target_link_libraries(time_simulation PRIVATE conway_lib)
target_link_libraries(time_simulation PRIVATE timing_lib)
target_link_libraries(time_count_neighbours PRIVATE conway_lib)
target_link_libraries(time_count_neighbours PRIVATE timing_lib)
target_link_libraries(parallel PRIVATE conway_lib)
target_link_libraries(parallel PRIVATE timing_lib)
target_link_libraries(time_1d_decomp PRIVATE conway_lib)
target_link_libraries(time_1d_decomp PRIVATE timing_lib)

find_package(MPI REQUIRED)
target_link_libraries(parallel PUBLIC MPI::MPI_CXX)
target_link_libraries(time_1d_decomp PUBLIC MPI::MPI_CXX)
