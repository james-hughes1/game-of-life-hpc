#include <iostream>
#include <string>

#include "conway/include/matrix.h"
#include "conway/include/world.h"

int main() {
    /**
     * Main procedure to run.
     */

    Matrix random_seed = matrix::generate_matrix(100, 100);
    World RandomWorld(random_seed);

    for (int age = 0; age < 100; age++) {
        RandomWorld.update_boundary();
        RandomWorld.evaluate_rules();
    }

    return 0;
}
