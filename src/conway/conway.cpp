/**
 * \file conway.cpp Defines the World class for simulating Game of Life.
 */

#include <random>

#include "conway.h"
#include "matrix.h"

World::World(int n_rows, int n_cols) {
    this->n_rows = n_rows;
    this->n_cols = n_cols;
    this->data   = new int[n_rows * n_cols * 2];
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, 1);
    for (int i = 0; i < n_rows * n_cols; i++) {
        this->data[i] = distrib(gen);
    }
    for (int i = n_rows * n_cols; i < n_rows * n_cols * 2; i++) {
        this->data[i] = 0;
    }
}

int &World::operator()(int i, int j) { return this->data[i * n_cols + j]; }

int World::evolve() { return 0; }

int World::update_boundary() { return 0; }

World::~World() { delete[] data; }

World::World(World &old_world) {
    n_rows = old_world.n_rows;
    n_cols = old_world.n_cols;
    data   = new int[n_rows * n_cols];
    for (int i = 0; i < n_rows * n_cols; i++)
        data[i] = old_world.data[i];
}

World &World::operator=(World &old_world) {
    n_rows = old_world.n_rows;
    n_cols = old_world.n_cols;
    data   = new int[n_rows * n_cols];
    for (int i = 0; i < n_rows * n_cols; i++)
        data[i] = old_world.data[i];
    return *this;
}
