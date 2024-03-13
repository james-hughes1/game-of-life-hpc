/**
 * \file main.cpp Assesses how the time taken for the simulation scales with
 * size.
 */

#include <fstream>
#include <iostream>
#include <string>

#include "conway/include/matrix.h"
#include "conway/include/world.h"
#include "timing/include/timing.h"

double time_simulation(int n_rows, int n_cols, int final_age) {
    Matrix random_seed = matrix::generate_matrix(n_rows, n_cols);
    World RandomWorld(random_seed);

    timing::start_clock();
    for (int age = 0; age < final_age; age++) {
        RandomWorld.update_boundary();
        RandomWorld.evaluate_rules();
    }
    std::cout << "Test entry: " << RandomWorld.Cells_0(1, 1) << std::endl;
    return timing::get_split();
}

int main() {
    /**
     * Main procedure to run.
     */

    std::fstream file("prof/time_simulation.txt");
    for (int world_size = 4; world_size <= 2048; world_size *= 2) {
        for (int offset = 0; offset < 3; offset++) {
            file << "World size " << world_size - offset << "x"
                 << world_size - offset << ", age 100 /////// Took: "
                 << time_simulation(world_size - offset, world_size - offset,
                                    100)
                 << " ms." << std::endl;
        }
        file << std::endl;
    }
    file.close();
    return 0;
}
