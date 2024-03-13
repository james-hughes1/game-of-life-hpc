/**
 * \file main.cpp Main program file.
 */

#include <iostream>
#include <string>

#include "matrix.h"
#include "world.h"

int main() {
    /**
     * Main procedure to run.
     */

    std::string glider_seed_str = matrix::read_file("seed_glider.txt");
    Matrix glider_seed          = matrix::read_matrix_str(glider_seed_str);
    World Glider(glider_seed);

    Glider.display_world();
    for (int age = 0; age < 50; age++) {
        Glider.update_boundary();
        Glider.evaluate_rules();
    }
    Glider.display_world();
    return 0;
}
