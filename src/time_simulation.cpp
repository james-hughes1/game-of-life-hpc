/**
 * \file main.cpp Assesses how the time taken for the simulation scales with
 * size.
 */

#include <iostream>
#include <string>

#include "conway/include/matrix.h"
#include "conway/include/world.h"
#include "timing/include/timing.h"

int time(int n_rows, int n_cols, int final_age) {
    Matrix random_seed = matrix::generate_matrix(n_rows, n_cols);
    World RandomWorld(random_seed);

    timing::start_clock();
    for (int age = 0; age < 100; age++) {
        RandomWorld.update_boundary();
        RandomWorld.evaluate_rules();
    }
    std::cout << "Test entry: " << RandomWorld.Cells_0(1, 1) << std::endl;
    std::cout << "World size " << n_rows << "x" << n_cols << " age "
              << final_age << " /////// Took: " << timing::get_split() << " ms."
              << std::endl;
    return 0;
}

int main() {
    /**
     * Main procedure to run.
     */

    for (int world_size = 2; world_size <= 1024; world_size *= 2) {
        time(world_size, world_size, 100);
    }
    for (int world_size = 4; world_size <= 1024; world_size *= 2) {
        time(world_size - 1, world_size - 1, 100);
    }
    for (int world_size = 4; world_size <= 1024; world_size *= 2) {
        time(world_size - 2, world_size - 2, 100);
    }

    return 0;
}
