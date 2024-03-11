/**
 * \file multiply.cpp Contains useful functions for matrix operations.
 */

#include <iostream>
#include <random>

#include "matrix.h"

Matrix::Matrix(int n_rows, int n_cols) {
    this->n_rows = n_rows;
    this->n_cols = n_cols;
    this->data   = new int[n_rows * n_cols];
}

int &Matrix::operator()(int i, int j) { return this->data[i * n_cols + j]; }

Matrix::~Matrix() { delete[] data; }

Matrix::Matrix(const Matrix &old_matrix) {
    n_rows = old_matrix.n_rows;
    n_cols = old_matrix.n_cols;
    data   = new int[n_rows * n_cols];
    for (int i = 0; i < n_rows * n_cols; i++) {
        data[i] = old_matrix.data[i];
    }
}
Matrix &Matrix::operator=(Matrix &old_matrix) {
    n_rows = old_matrix.n_rows;
    n_cols = old_matrix.n_cols;
    data   = new int[n_rows * n_cols];
    for (int i = 0; i < n_rows * n_cols; i++) {
        data[i] = old_matrix.data[i];
    }
    return *this;
}

int matrix::display_matrix(Matrix A) {
    /**
     * \brief Displays a matrix.
     * \param A Matrix to be displayed.
     */
    for (int i = 0; i < A.n_rows; i++) {
        for (int j = 0; j < A.n_cols; j++) {
            std::cout << A(i, j) << " ";
        }
        std::cout << std::endl;
    }
    return 0;
}

Matrix matrix::count_neighbours(Matrix A) {
    /**
     * \brief Performs the necessary 3x3 2D convolution required for Game of
     * Life. \param A Input matrix.
     */

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
    return B;
}

Matrix matrix::generate_matrix(int n_rows, int n_cols) {
    Matrix A(n_rows, n_cols);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, 1);
    for (int i = 0; i < n_rows; i++) {
        for (int j = 0; j < n_cols; j++) {
            A(i, j) = distrib(gen);
        }
    }
    return A;
}
