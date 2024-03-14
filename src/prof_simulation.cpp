#include <iostream>
#include <string>

#include "conway/include/matrix.h"
#include "conway/include/world.h"

int main() {
    /**
     * Main procedure to run.
     */

    Matrix random_seed = matrix::generate_matrix(1000, 1000);
    World RandomWorld(random_seed);

    for (int age = 0; age < 10; age++) {
        RandomWorld.update_boundary();
        RandomWorld.evaluate_rules();
    }

    std::cout << "Test entry: " << RandomWorld.Cells_0(1, 1) << std::endl;

    return 0;
}
