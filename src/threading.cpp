#include <iostream>
#include <omp.h>

#include "conway/include/matrix.h"

int main() {

    int N_ROWS   = 8000;
    int N_COLS   = 8000;
    Matrix cells = matrix::generate_matrix(N_ROWS, N_COLS);
    Matrix row_convolution(N_ROWS, N_COLS), counts(N_ROWS, N_COLS);

#pragma omp parallel for
    // Row convolution
    for (int i = 0; i < N_ROWS; i++) {
        for (int j = 1; j < N_COLS - 1; j++) {
            row_convolution(i, j) =
                cells(i, j - 1) + cells(i, j) + cells(i, j + 1);
        }
    }

#pragma omp barrier

#pragma omp parallel for
    // Column convolution
    for (int i = 1; i < N_ROWS - 1; i++) {
        for (int j = 1; j < N_COLS - 1; j++) {
            counts(i, j) = row_convolution(i - 1, j) + row_convolution(i, j) +
                           row_convolution(i + 1, j);
        }
    }

    std::cout << "Test entry: " << counts(1, 1) << std::endl;
}
