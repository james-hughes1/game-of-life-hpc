add_subdirectory(conway)

add_executable(main main.cpp)
ClangFormat(main .)

target_link_libraries(main PRIVATE conway_lib)
