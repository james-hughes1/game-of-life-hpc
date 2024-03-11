/**
 * \file conway.cpp Defines the World class for simulating Game of Life.
 */

#include <random>

#include "matrix.h"
#include "world.h"

World::World(Matrix seed)
    : cells_0(seed.n_rows, seed.n_cols), cells_1(seed.n_rows, seed.n_cols) {
    this->cells_0 = seed;
    this->cells_1 = seed;
    this->n_rows  = cells_0.n_rows;
    this->n_cols  = cells_0.n_cols;
    this->age     = 0;
    this->side    = 0;
}

int World::evolve() { return 0; }

int World::update_boundary() { return 0; }

int World::display_world() {
    if (this->side == 0) {
        matrix::display_matrix(this->cells_0);
    } else {
        matrix::display_matrix(this->cells_1);
    }
    return 0;
}
