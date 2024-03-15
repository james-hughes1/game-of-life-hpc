#include <iostream>
#include <omp.h>

#include "matrix.h"

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
#pragma omp parallel for
    // Row convolution
    for (int i = 0; i < n_rows; i++) {
        for (int j = 1; j < n_cols - 1; j++) {
            row_convolution(i, j) =
                cells_0(i, j - 1) + cells_0(i, j) + cells_0(i, j + 1);
        }
    }

#pragma omp barrier

#pragma omp parallel for
    // Column convolution
    for (int i = 1; i < n_rows - 1; i++) {
        for (int j = 1; j < n_cols - 1; j++) {
            counts(i, j) = row_convolution(i - 1, j) + row_convolution(i, j) +
                           row_convolution(i + 1, j);
        }
    }

#pragma omp barrier

#pragma omp parallel for
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
    omp_set_num_threads(4);

    int N_ROWS           = 43;
    int N_COLS           = 17;
    int MAX_AGE          = 50;
    std::string seed_str = matrix::read_file("test/test_data/input_file_2.txt");
    Matrix seed          = matrix::read_matrix_str(seed_str);
    std::cout << "Initial seed:\n" << matrix::write_matrix_str(seed);

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

    Matrix final_state(N_ROWS, N_COLS);
    if (age % 2 == 0) {
        for (int i = 0; i < N_ROWS; i++) {
            for (int j = 0; j < N_COLS; j++) {
                final_state(i, j) = cells_0(i + 1, j + 1);
            }
        }
    } else {
        for (int i = 0; i < N_ROWS; i++) {
            for (int j = 0; j < N_COLS; j++) {
                final_state(i, j) = cells_1(i + 1, j + 1);
            }
        }
    }
    std::cout << "Final state at age " << MAX_AGE << ":\n"
              << matrix::write_matrix_str(final_state);
    std::string expected_str =
        matrix::read_file("test/test_data/output_file_2.txt");
    Matrix expected_state = matrix::read_matrix_str(expected_str);
    if (final_state == expected_state) {
        std::cout << "Test passed." << std::endl;
    } else {
        std::cout << "Test failed." << std::endl;
    }

    return 0;
}
