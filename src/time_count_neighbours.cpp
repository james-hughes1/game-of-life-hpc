/**
 * \file main.cpp Assesses how the time taken for the convolution scales with
 * size.
 */

#include <fstream>
#include <iostream>
#include <string>

#include "conway/include/matrix.h"
#include "timing/include/timing.h"

double time_count_neighbours(int n_rows, int n_cols) {
    double total_time = 0;
    for (int run = 0; run < 10; run++) {
        Matrix random_matrix = matrix::generate_matrix(n_rows, n_cols);
        timing::start_clock();
        Matrix count_matrix = matrix::count_neighbours(random_matrix);
        total_time += timing::get_split();
    }
    return total_time / 10;
}

int main() {
    /**
     * Main procedure to run.
     */

    std::fstream file("prof/time_count_neighbours.txt");
    for (int matrix_size = 4; matrix_size <= 10000; matrix_size *= 2) {
        for (int offset = 0; offset < 3; offset++) {
            file << "Matrix size " << matrix_size - offset << "x"
                 << matrix_size - offset << " /////// Took: "
                 << time_count_neighbours(matrix_size - offset,
                                          matrix_size - offset)
                 << " ms. (Mean of 10 runs.)" << std::endl;
        }
        file << std::endl;
    }
    file.close();

    return 0;
}
