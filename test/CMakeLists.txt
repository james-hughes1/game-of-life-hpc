add_executable(test_conway
    test_matrix.cpp
    test_world.cpp
)
ClangFormat(test_conway .)

add_executable(test_dd1
    test_dd1.cpp
)
ClangFormat(test_dd1 .)

add_executable(test_dd2
    test_dd2.cpp
)
ClangFormat(test_dd2 .)

add_executable(test_omp
    test_omp.cpp
)
ClangFormat(test_omp .)

add_executable(test_hybrid
    test_hybrid.cpp
)
ClangFormat(test_hybrid .)

target_link_libraries(test_conway
    PUBLIC
    conway_lib
    GTest::gtest_main
)
target_link_libraries(test_dd1
    PUBLIC
    conway_lib
)
target_link_libraries(test_dd2
    PUBLIC
    conway_lib
)
target_link_libraries(test_omp
    PUBLIC
    conway_lib
)
target_link_libraries(test_hybrid
    PUBLIC
    conway_lib
)

find_package(MPI REQUIRED)
target_link_libraries(test_dd1 PUBLIC MPI::MPI_CXX)
target_link_libraries(test_dd2 PUBLIC MPI::MPI_CXX)
target_link_libraries(test_hybrid PUBLIC MPI::MPI_CXX)

find_package(OpenMP)
target_link_libraries(test_omp PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(test_hybrid PUBLIC OpenMP::OpenMP_CXX)

target_compile_definitions(test_omp PRIVATE fopenmp O3 march=native)
target_compile_definitions(test_hybrid PRIVATE fopenmp O3 march=native)

include(GoogleTest)
gtest_discover_tests(test_conway)
