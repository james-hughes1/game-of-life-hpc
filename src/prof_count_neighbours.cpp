#include <iostream>
#include <string>

#include "conway/include/matrix.h"

int main() {
    /**
     * \brief Performs the necessary 3x3 2D convolution required for Game of
     * Life. \param A Input matrix.
     */

    Matrix A = matrix::generate_matrix(1000, 1000);

    int n_rows = A.n_rows;
    int n_cols = A.n_cols;

    // Row convolution
    Matrix B(n_rows, n_cols);
    for (int i = 0; i < n_rows; i++) {
        for (int j = 1; j < n_cols - 1; j++) {
            B(i, j) = A(i, j - 1) + A(i, j) + A(i, j + 1);
        }
    }

    // Column convolution (and subtract 1 for the proper convolution value).
    Matrix C(n_rows, n_cols);
    for (int i = 1; i < n_rows - 1; i++) {
        for (int j = 1; j < n_cols - 1; j++) {
            C(i, j) = B(i - 1, j) + B(i, j) + B(i + 1, j) - A(i, j);
        }
    }

    std::cout << "Test entry: " << C(1, 1) << std::endl;

    return 0;
}
