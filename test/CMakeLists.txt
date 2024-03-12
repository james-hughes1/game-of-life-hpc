add_executable(test_conway
    test_matrix.cpp
    test_world.cpp
)
ClangFormat(test_conway .)

target_link_libraries(test_conway
    PUBLIC
    conway_lib
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(test_conway)
