/**
 * \file conway.cpp Defines the World class for simulating Game of Life.
 */

#include <iostream>
#include <random>
#include <string>

#include "matrix.h"
#include "world.h"

World::World(Matrix seed) : Cells_0(seed), Cells_1(seed.n_rows, seed.n_cols) {
    n_rows = seed.n_rows;
    n_cols = seed.n_cols;
    age    = 0;
}

int World::evaluate_rules() {
    if (age % 2 == 0) {
        Matrix Cells_count = matrix::count_neighbours(Cells_0);
        conway::evaluate_rules(Cells_count, Cells_0, Cells_1);
    } else {
        Matrix Cells_count = matrix::count_neighbours(Cells_1);
        conway::evaluate_rules(Cells_count, Cells_1, Cells_0);
    }
    age += 1;
    return 0;
}

int World::update_boundary() {
    // Vertices
    if (age % 2 == 0) {
        conway::update_boundary(Cells_0);
    } else {
        conway::update_boundary(Cells_1);
    }
    return 0;
}

int World::display_world() {
    std::cout << "World has age: " << age << " ticks." << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::string output = (age % 2 == 0) ? matrix::write_matrix(Cells_0)
                                        : matrix::write_matrix(Cells_1);
    std::cout << output << std::endl;
    return 0;
}

Matrix World::output_cells() { return Cells_0; }

int conway::evaluate_rules(Matrix Cells_count, Matrix &Cells_current,
                           Matrix &Cells_next) {
    for (int i = 1; i < Cells_count.n_rows - 1; i++) {
        for (int j = 1; j < Cells_count.n_cols - 1; j++) {
            if (Cells_current(i, j) == 1) {
                if (Cells_count(i, j) != 2 and Cells_count(i, j) != 3) {
                    Cells_next(i, j) = 0;
                } else {
                    Cells_next(i, j) = 1;
                }
            } else {
                if (Cells_count(i, j) == 3) {
                    Cells_next(i, j) = 1;
                } else {
                    Cells_next(i, j) = 0;
                }
            }
        }
    }
    return 0;
}

int conway::update_boundary(Matrix &Cells_current) {
    // Vertices
    Cells_current(Cells_current.n_rows - 1, Cells_current.n_cols - 1) =
        Cells_current(1, 1);
    Cells_current(0, Cells_current.n_cols - 1) =
        Cells_current(Cells_current.n_rows - 2, 1);
    Cells_current(0, 0) =
        Cells_current(Cells_current.n_rows - 2, Cells_current.n_cols - 2);
    Cells_current(Cells_current.n_rows - 1, 0) =
        Cells_current(1, Cells_current.n_cols - 2);
    // Edges
    for (int j = 1; j < Cells_current.n_cols - 1; j++) {
        Cells_current(0, j) = Cells_current(Cells_current.n_rows - 2, j);
        Cells_current(Cells_current.n_rows - 1, j) = Cells_current(1, j);
    }
    for (int i = 1; i < Cells_current.n_cols - 1; i++) {
        Cells_current(i, 0) = Cells_current(i, Cells_current.n_cols - 2);
        Cells_current(i, Cells_current.n_rows - 1) = Cells_current(i, 1);
    }
    return 0;
}
