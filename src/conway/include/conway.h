#ifndef CONWAY_H
#define CONWAY_H

#include "matrix.h"

/**
 * @brief Library for the functions involved in simulating a Game of Life world.
 */

class World {
  private:
    int *data;

  public:
    int n_rows;
    int n_cols;
    int age;
    int side;
    World(int n_rows, int n_cols);
    int &operator()(int i, int j);
    int evolve();
    int update_boundary();
    ~World();
    World(World &old_world);
    World &operator=(World &world);
};

#endif
