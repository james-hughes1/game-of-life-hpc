#include <iostream>
#include <string>

#include "conway/include/matrix.h"

int main() {
    /**
     * \brief Performs the necessary 3x3 2D convolution required for Game of
     * Life. \param A Input matrix.
     */

    Matrix A = matrix::generate_matrix(1000, 1000);

    Matrix B(A.n_rows, A.n_cols);

    for (int i = 0; i < A.n_rows; i++) {
        for (int j = 0; j < A.n_cols; j++) {
            B(i, j) = 0;
            for (int p = -1; p < 2; p++) {
                for (int q = -1; q < 2; q++) {
                    B(i, j) += A(i + p, j + q);
                }
            }
            B(i, j) -= A(i, j);
            if (i == 0 or i == A.n_rows - 1 or j == 0 or j == A.n_cols - 1) {
                B(i, j) = -1;
            }
        }
    }

    return 0;
}
