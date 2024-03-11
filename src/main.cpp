/**
 * \file main.cpp Main program file.
 */

#include <iostream>

#include "matrix.h"

int main() {
    /**
     * Main procedure to run.
     */
    Matrix A = matrix::generate_matrix(4, 4);
    matrix::display_matrix(A);
    Matrix B = matrix::count_neighbours(A);
    matrix::display_matrix(B);
    return 0;
}
