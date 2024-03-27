/**
 * \file main.cpp Main program file.
 */

#include <iostream>
#include <string>

#include "conway/include/matrix.h"
#include "conway/include/world.h"

int main(int argc, char **argv) {
    /**
     * Main procedure to run.
     */

    int max_age = std::stoi(argv[1]);
    Matrix seed =
        (argc == 4)
            ? matrix::generate_matrix(std::stoi(argv[2]), std::stoi(argv[3]))
            : matrix::read_matrix_str(matrix::read_file(argv[2]));

    World InputWorld(seed);

    InputWorld.display_world();
    for (int age = 0; age < max_age; age++) {
        InputWorld.update_boundary();
        InputWorld.evaluate_rules();
    }
    InputWorld.display_world();
    return 0;
}
