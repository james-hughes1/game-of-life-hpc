/**
 * \file conway.cpp Defines the World class for simulating Game of Life.
 */

#include <iostream>
#include <random>

#include "matrix.h"
#include "world.h"

World::World(Matrix seed) : cells_0(seed), cells_1(seed.n_rows, seed.n_cols) {
    n_rows = seed.n_rows;
    n_cols = seed.n_cols;
    age    = 0;
}

int World::evaluate_rules() {
    Matrix cells_count = (age % 2 == 0) ? matrix::count_neighbours(cells_0)
                                        : matrix::count_neighbours(cells_1);
    matrix::display_matrix(cells_count);
    Matrix *cells_old;
    Matrix *cells_new;
    cells_old = (age % 2 == 0) ? &cells_0 : &cells_1;
    cells_new = (age % 2 == 0) ? &cells_1 : &cells_0;
    conway::evaluate_rules(cells_count, cells_old, cells_new);
    age += 1;
    return 0;
}

int World::update_boundary() { return 0; }

int World::display_world() {
    std::cout << "World has age: " << age << " ticks..." << std::endl;
    if (age % 2 == 0) {
        matrix::display_matrix(cells_0);
    } else {
        matrix::display_matrix(cells_1);
    }
    return 0;
}

int conway::evaluate_rules(Matrix cells_count, Matrix *cells_old,
                           Matrix *cells_new) {
    for (int i = 1; i < cells_count.n_rows - 1; i++) {
        for (int j = 1; j < cells_count.n_cols - 1; j++) {
            if ((*cells_old)(i, j) == 1) {
                if (cells_count(i, j) != 2 and cells_count(i, j) != 3) {
                    (*cells_new)(i, j) = 0;
                } else {
                    (*cells_new)(i, j) = 1;
                }
            } else {
                if (cells_count(i, j) == 3) {
                    (*cells_new)(i, j) = 1;
                } else {
                    (*cells_new)(i, j) = 0;
                }
            }
        }
    }
    return 0;
}
