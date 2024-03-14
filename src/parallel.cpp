/**
 * \file main.cpp Main program file.
 */

#include <iostream>
#include <string>
#include <tuple>

#include "conway/include/matrix.h"
#include "conway/include/parallel.h"
#include "conway/include/world.h"

int main(int argc, char **argv) {
    /**
     * Main procedure to run.
     */

    parallel::decomposition_1d(argc, argv);
}
