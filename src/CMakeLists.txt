add_subdirectory(conway)
add_subdirectory(timing)

add_executable(main main.cpp)
ClangFormat(main .)

add_executable(time_simulation time_simulation.cpp)
ClangFormat(time_simulation .)

target_link_libraries(main PRIVATE conway_lib)
target_link_libraries(time_simulation PRIVATE conway_lib)
target_link_libraries(time_simulation PRIVATE timing_lib)
