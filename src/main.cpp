/**
 * \file main.cpp Main program file.
 */

#include <iostream>

#include "matrix.h"
#include "world.h"

int main() {
    /**
     * Main procedure to run.
     */

    Matrix A = matrix::generate_matrix(6, 10);

    World world(A);

    world.display_world();
    world.update_boundary();
    world.display_world();
    world.evaluate_rules();
    world.display_world();
    world.update_boundary();
    world.display_world();

    matrix::read_file("test/test_data/input_file_1.txt");

    return 0;
}
