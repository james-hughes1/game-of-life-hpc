add_executable(test_conway
    test_matrix.cpp
    test_world.cpp
)
ClangFormat(test_conway .)
add_executable(test_parallel
    test_parallel.cpp
)
ClangFormat(test_parallel .)

target_link_libraries(test_conway
    PUBLIC
    conway_lib
    GTest::gtest_main
)
target_link_libraries(test_parallel
    PUBLIC
    conway_lib
)

find_package(MPI REQUIRED)
target_link_libraries(test_parallel PUBLIC MPI::MPI_CXX)

include(GoogleTest)
gtest_discover_tests(test_conway)
