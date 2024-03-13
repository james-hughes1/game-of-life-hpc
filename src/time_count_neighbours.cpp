/**
 * \file main.cpp Assesses how the time taken for the convolution scales with
 * size.
 */

#include <iostream>
#include <string>

#include "conway/include/matrix.h"
#include "timing/include/timing.h"

int time_count_neighbours(int n_rows, int n_cols) {
    double total_time;
    timing::start_clock();
    for (int run = 0; run < 10; run++) {
        Matrix random_matrix = matrix::generate_matrix(n_rows, n_cols);
        double start_time    = timing::get_split();
        matrix::count_neighbours(random_matrix);
        total_time += (timing::get_split() - start_time);
    }
    std::cout << "Matrix size " << n_rows << "x" << n_cols
              << " /////// Took: " << timing::get_split() / 10
              << " ms. (Mean of 10 runs.)" << std::endl;
    return 0;
}

int main() {
    /**
     * Main procedure to run.
     */

    for (int matrix_size = 2; matrix_size <= 1024; matrix_size *= 2) {
        time_count_neighbours(matrix_size, matrix_size);
    }
    for (int matrix_size = 4; matrix_size <= 1024; matrix_size *= 2) {
        time_count_neighbours(matrix_size - 1, matrix_size - 1);
    }
    for (int matrix_size = 4; matrix_size <= 1024; matrix_size *= 2) {
        time_count_neighbours(matrix_size - 2, matrix_size - 2);
    }

    return 0;
}
