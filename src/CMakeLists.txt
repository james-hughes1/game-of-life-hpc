add_subdirectory(conway)
add_subdirectory(timing)

add_executable(run_single run_single.cpp)
ClangFormat(run_single .)

add_executable(time_single time_single.cpp)
ClangFormat(time_single .)

add_executable(time_count_neighbours time_count_neighbours.cpp)
ClangFormat(time_count_neighbours .)

add_executable(run_dd1 run_dd1.cpp)
ClangFormat(run_dd1 .)

add_executable(time_dd1
time_dd1.cpp
)
ClangFormat(time_dd1 .)

add_executable(run_dd2
run_dd2.cpp
)
ClangFormat(run_dd2 .)

add_executable(time_dd2
time_dd2.cpp
)
ClangFormat(time_dd2 .)

add_executable(run_omp
run_omp.cpp
)
ClangFormat(run_omp .)

add_executable(run_hybrid
run_hybrid.cpp
)
ClangFormat(run_hybrid .)

add_executable(time_hybrid
time_hybrid.cpp
)
ClangFormat(time_hybrid .)

target_link_libraries(run_single PRIVATE conway_lib)
target_link_libraries(time_single PRIVATE conway_lib)
target_link_libraries(time_single PRIVATE timing_lib)
target_link_libraries(time_count_neighbours PRIVATE conway_lib)
target_link_libraries(time_count_neighbours PRIVATE timing_lib)
target_link_libraries(run_dd1 PRIVATE conway_lib)
target_link_libraries(run_dd1 PRIVATE timing_lib)
target_link_libraries(time_dd1 PRIVATE conway_lib)
target_link_libraries(time_dd1 PRIVATE timing_lib)
target_link_libraries(run_dd2 PRIVATE conway_lib)
target_link_libraries(run_dd2 PRIVATE timing_lib)
target_link_libraries(time_dd2 PRIVATE conway_lib)
target_link_libraries(time_dd2 PRIVATE timing_lib)
target_link_libraries(run_omp PRIVATE conway_lib)
target_link_libraries(run_omp PRIVATE timing_lib)
target_link_libraries(run_hybrid PRIVATE conway_lib)
target_link_libraries(run_hybrid PRIVATE timing_lib)
target_link_libraries(time_hybrid PRIVATE conway_lib)
target_link_libraries(time_hybrid PRIVATE timing_lib)

find_package(MPI REQUIRED)
target_link_libraries(run_dd1 PUBLIC MPI::MPI_CXX)
target_link_libraries(time_dd1 PUBLIC MPI::MPI_CXX)
target_link_libraries(run_dd2 PUBLIC MPI::MPI_CXX)
target_link_libraries(time_dd2 PUBLIC MPI::MPI_CXX)
target_link_libraries(run_hybrid PUBLIC MPI::MPI_CXX)
target_link_libraries(time_hybrid PUBLIC MPI::MPI_CXX)

find_package(OpenMP)
target_link_libraries(run_omp PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(run_hybrid PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(time_hybrid PUBLIC OpenMP::OpenMP_CXX)

target_compile_definitions(run_omp PRIVATE fopenmp O3 march=native)
target_compile_definitions(run_hybrid PRIVATE fopenmp O3 march=native)
target_compile_definitions(time_hybrid PRIVATE fopenmp O3 march=native)
