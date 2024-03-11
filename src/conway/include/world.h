#ifndef CONWAY_H
#define CONWAY_H

#include "matrix.h"

/**
 * @brief Library for the functions involved in simulating a Game of Life world.
 */

class World {
  public:
    int n_rows;
    int n_cols;
    Matrix cells_0;
    Matrix cells_1;
    int age;
    int side;
    World(Matrix seed);
    int evolve();
    int update_boundary();
    int display_world();
    int random_seed();
};

#endif
