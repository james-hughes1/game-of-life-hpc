#include <iostream>
#include <omp.h>

#include "conway/include/matrix.h"

int update_boundary_omp(Matrix &cells_0) {
    // Update vertices
    cells_0(cells_0.n_rows - 1, cells_0.n_cols - 1) = cells_0(1, 1);
    cells_0(0, cells_0.n_cols - 1) = cells_0(cells_0.n_rows - 2, 1);
    cells_0(0, 0) = cells_0(cells_0.n_rows - 2, cells_0.n_cols - 2);
    cells_0(cells_0.n_rows - 1, 0) = cells_0(1, cells_0.n_cols - 2);
    // Update Edges

#pragma omp parallel for
    for (int j = 1; j < cells_0.n_cols - 1; j++) {
        cells_0(0, j)                  = cells_0(cells_0.n_rows - 2, j);
        cells_0(cells_0.n_rows - 1, j) = cells_0(1, j);
    }

#pragma omp parallel for
    for (int i = 1; i < cells_0.n_rows - 1; i++) {
        cells_0(i, 0)                  = cells_0(i, cells_0.n_cols - 2);
        cells_0(i, cells_0.n_cols - 1) = cells_0(i, 1);
    }
    return 0;
}

int evolve_omp(Matrix &cells_0, Matrix &cells_1) {
    int n_rows = cells_0.n_rows;
    int n_cols = cells_0.n_cols;
    Matrix row_convolution(n_rows, n_cols);
    Matrix counts(n_rows, n_cols);
#pragma omp parallel for collapse(2)
    // Row convolution
    for (int i = 0; i < n_rows; i++) {
        for (int j = 1; j < n_cols - 1; j++) {
            row_convolution(i, j) =
                cells_0(i, j - 1) + cells_0(i, j) + cells_0(i, j + 1);
        }
    }

#pragma omp barrier

#pragma omp parallel for collapse(2)
    // Column convolution
    for (int i = 1; i < n_rows - 1; i++) {
        for (int j = 1; j < n_cols - 1; j++) {
            counts(i, j) = row_convolution(i - 1, j) + row_convolution(i, j) +
                           row_convolution(i + 1, j);
        }
    }

#pragma omp barrier

#pragma omp parallel for collapse(2)
    // Evaluate rules
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            cells_1(i, j) = (counts(i, j) == 3) ||
                            ((counts(i, j) == 4) && (cells_0(i, j) == 1));
        }
    }
    return 0;
}

int main() {

    int N_ROWS  = 4000;
    int N_COLS  = 4000;
    int MAX_AGE = 50;
    Matrix seed = matrix::generate_matrix(N_ROWS, N_COLS);

    Matrix cells_0(N_ROWS + 2, N_COLS + 2);
    Matrix cells_1(N_ROWS + 2, N_COLS + 2);
    cells_0.write_sub_matrix(seed);

    int age = 0;
    while (age < MAX_AGE) {
        // Update and evaluate rules.
        if (age % 2 == 0) {
            update_boundary_omp(cells_0);
            evolve_omp(cells_0, cells_1);
        } else {
            update_boundary_omp(cells_1);
            evolve_omp(cells_1, cells_0);
        }
        age++;
    }

    if (age % 2 == 0) {
        std::cout << "test " << cells_0(1, 1) << std::endl;
    } else {
        std::cout << "test " << cells_1(1, 1) << std::endl;
    }

    return 0;
}
