/**
 * \file multiply.cpp Contains useful functions for matrix operations.
 */

#include <iostream>
#include <random>
#include <string>

#include "matrix.h"

// class Matrix

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

bool Matrix::operator==(Matrix &matrix_other) {
    if (this->n_rows != matrix_other.n_rows or
        this->n_cols != matrix_other.n_cols) {
        return false;
    }
    for (int i = 0; i < (this->n_rows) * (this->n_cols); i++) {
        if (this->data[i] != matrix_other.data[i]) {
            return false;
        }
    }
    return true;
}

int Matrix::zero() {
    for (int i = 0; i < (n_rows) * (n_cols); i++) {
        data[i] = 0;
    }
    return 0;
}

int Matrix::write_submatrix(Matrix &submatrix) {
    this->data[1] = submatrix(0, 0);
    return 0;
}
Matrix Matrix::read_submatrix() { return *this; }

// namespace matrix

std::string matrix::read_file(std::string filename) { return filename; }

Matrix matrix::read_matrix(std::string filename) {
    std::cout << filename[0] << std::endl;
    Matrix A(1, 1);
    return A;
}

std::string matrix::write_matrix(Matrix A) {
    std::cout << A(0, 0) << std::endl;
    std::string output_str = "X";
    return output_str;
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
