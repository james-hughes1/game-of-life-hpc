add_executable(test_conway
    test_matrix.cpp
    test_world.cpp
)
ClangFormat(test_conway .)

add_executable(test_parallel
    test_parallel.cpp
)
ClangFormat(test_parallel .)

add_executable(test_parallel_2d
    test_parallel_2d.cpp
)
ClangFormat(test_parallel_2d .)

add_executable(test_threading
    test_threading.cpp
)
ClangFormat(test_threading .)

target_link_libraries(test_conway
    PUBLIC
    conway_lib
    GTest::gtest_main
)
target_link_libraries(test_parallel
    PUBLIC
    conway_lib
)
target_link_libraries(test_parallel_2d
    PUBLIC
    conway_lib
)
target_link_libraries(test_threading
    PUBLIC
    conway_lib
)

find_package(MPI REQUIRED)
target_link_libraries(test_parallel PUBLIC MPI::MPI_CXX)
target_link_libraries(test_parallel_2d PUBLIC MPI::MPI_CXX)

find_package(OpenMP)
target_link_libraries(test_threading PUBLIC OpenMP::OpenMP_CXX)

target_compile_definitions(test_threading PRIVATE fopenmp O3 march=native)

include(GoogleTest)
gtest_discover_tests(test_conway)
