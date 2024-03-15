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

add_executable(parallel_2d
parallel_2d.cpp
)
ClangFormat(parallel_2d .)

add_executable(time_2d_decomp
time_2d_decomp.cpp
)
ClangFormat(time_2d_decomp .)

add_executable(threading
threading.cpp
)
ClangFormat(threading .)

target_link_libraries(main PRIVATE conway_lib)
target_link_libraries(time_simulation PRIVATE conway_lib)
target_link_libraries(time_simulation PRIVATE timing_lib)
target_link_libraries(time_count_neighbours PRIVATE conway_lib)
target_link_libraries(time_count_neighbours PRIVATE timing_lib)
target_link_libraries(parallel PRIVATE conway_lib)
target_link_libraries(parallel PRIVATE timing_lib)
target_link_libraries(time_1d_decomp PRIVATE conway_lib)
target_link_libraries(time_1d_decomp PRIVATE timing_lib)
target_link_libraries(parallel_2d PRIVATE conway_lib)
target_link_libraries(parallel_2d PRIVATE timing_lib)
target_link_libraries(time_2d_decomp PRIVATE conway_lib)
target_link_libraries(time_2d_decomp PRIVATE timing_lib)
target_link_libraries(threading PRIVATE conway_lib)
target_link_libraries(threading PRIVATE timing_lib)

find_package(MPI REQUIRED)
target_link_libraries(parallel PUBLIC MPI::MPI_CXX)
target_link_libraries(time_1d_decomp PUBLIC MPI::MPI_CXX)
target_link_libraries(parallel_2d PUBLIC MPI::MPI_CXX)
target_link_libraries(time_2d_decomp PUBLIC MPI::MPI_CXX)

find_package(OpenMP)
target_link_libraries(threading PUBLIC OpenMP::OpenMP_CXX)

target_compile_definitions(threading PRIVATE fopenmp O3 march=native)
