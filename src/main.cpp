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

    Matrix A = matrix::generate_matrix(5, 5);

    World world(A);

    world.display_world();

    world.evaluate_rules();

    world.display_world();

    world.evaluate_rules();

    world.display_world();

    world.evaluate_rules();

    world.display_world();

    return 0;
}
