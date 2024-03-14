/**
 * \file conway.cpp Defines the World class for simulating Game of Life.
 */

#include <iostream>
#include <random>
#include <string>

#include "include/matrix.h"
#include "include/world.h"

World::World(Matrix seed)
    : Cells_0(seed.n_rows + 2, seed.n_cols + 2),
      Cells_1(seed.n_rows + 2, seed.n_cols + 2) {
    Cells_0.write_sub_matrix(seed);
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
    if (age % 2 == 0) {
        conway::update_boundary(Cells_0);
    } else {
        conway::update_boundary(Cells_1);
    }
    return 0;
}

Matrix World::output_cells() {
    if (age % 2 == 0) {
        return Cells_0.read_sub_matrix();
    } else {
        return Cells_1.read_sub_matrix();
    }
}

int World::display_world() {
    std::cout << "World has age: " << age << " ticks." << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::string output = matrix::write_matrix_str((*this).output_cells());
    std::cout << output << std::endl;
    return 0;
}

int World::write_edge(int *edge, int loc) {
    // loc==0 means write top, loc==1 means bottom.
    if (age % 2 == 0) {
        for (int j = 0; j < Cells_0.n_cols; j++) {
            Cells_0(loc * (Cells_0.n_rows - 1), j) = edge[j];
        }
    } else {
        for (int j = 0; j < Cells_1.n_cols; j++) {
            Cells_1(loc * (Cells_1.n_rows - 1), j) = edge[j];
        }
    }
    return 0;
}

int World::read_edge(int *edge, int loc) {
    // loc==0 means write top, loc==1 means bottom.
    if (age % 2 == 0) {
        for (int j = 0; j < Cells_0.n_cols; j++) {
            edge[j] = Cells_0(loc * (Cells_0.n_rows - 1), j);
        }
    } else {
        for (int j = 0; j < Cells_1.n_cols; j++) {
            edge[j] = Cells_1(loc * (Cells_1.n_rows - 1), j);
        }
    }
    return 0;
}

int conway::evaluate_rules(Matrix &Cells_count, Matrix &Cells_current,
                           Matrix &Cells_next) {
    for (int i = 0; i < Cells_count.n_rows; i++) {
        for (int j = 0; j < Cells_count.n_cols; j++) {
            Cells_next(i, j) =
                (Cells_count(i, j) == 3) ||
                ((Cells_count(i, j) == 4) && (Cells_current(i, j) == 1));
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
    for (int i = 1; i < Cells_current.n_rows - 1; i++) {
        Cells_current(i, 0) = Cells_current(i, Cells_current.n_cols - 2);
        Cells_current(i, Cells_current.n_cols - 1) = Cells_current(i, 1);
    }
    return 0;
}
