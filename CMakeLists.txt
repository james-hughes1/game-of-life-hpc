cmake_minimum_required(VERSION 3.15...3.26)

project(
    game_of_life
    VERSION 1.0
    LANGUAGES CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

list(APPEND CMAKE_MODULE_PATH "./cmake")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

set(CMAKE_CXX_FLAGS "-O3 -march=native -Werror -Wall -Wextra -pedantic")

include(Doxygen)
include(ClangFormat)
include(Testing)

set(DOXYGEN_PROJECT_NAME "Game of Life")
doxygen_add_docs(
    doxygen
    ${PROJECT_SOURCE_DIR}
    COMMENT "Generating documentation"
)

add_subdirectory(src)
add_subdirectory(test)
