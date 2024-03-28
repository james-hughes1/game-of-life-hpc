/**
 * @file time_single.cpp Assesses how the time taken for the simulation scales
 * with size, for 200 ticks.
 * @details Can be run by `./bin/time_single`
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
    return timing::get_split();
}

int main() {
    /**
     * Main procedure to run.
     */

    int MAX_AGE = 200;
    std::fstream file("prof/time_simulation.txt");
    for (int world_size = 1000; world_size <= 10000; world_size += 1000) {
        file << world_size << " " << world_size << " " << MAX_AGE << " "
             << time_simulation(world_size, world_size, MAX_AGE) << std::endl;
    }
    file.close();
    return 0;
}
