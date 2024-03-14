#ifndef CONWAY_H
#define CONWAY_H

#include <tuple>

#include "matrix.h"

/**
 * @brief Library for the functions involved in simulating a Game of Life world.
 */

class World {
  public:
    int n_rows;
    int n_cols;
    Matrix Cells_0;
    Matrix Cells_1;
    int age;
    int side;
    World(Matrix seed);
    int evaluate_rules();
    int update_boundary();
    int display_world();
    int random_seed();
    Matrix output_cells();
    int write_edge(int *edge, int loc);
    int read_edge(int *edge, int loc);
};

namespace conway {
int evaluate_rules(Matrix &Cells_count, Matrix &Cells_current,
                   Matrix &Cells_next);
int update_boundary(Matrix &Cells);
std::tuple<int, int> divide_rows(int rows, int size, int rank);
}; // namespace conway

#endif
