/**
 * \file multiply.cpp Contains useful functions for matrix operations.
 */

#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "include/matrix.h"

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

int Matrix::write_sub_matrix(Matrix &sub_matrix) {
    for (int i = 0; i < sub_matrix.n_rows; i++) {
        for (int j = 0; j < sub_matrix.n_cols; j++) {
            (*this)(i + 1, j + 1) = sub_matrix(i, j);
        }
    }
    return 0;
}

Matrix Matrix::read_sub_matrix() {
    Matrix sub_matrix(this->n_rows - 2, this->n_cols - 2);
    for (int i = 0; i < sub_matrix.n_rows; i++) {
        for (int j = 0; j < sub_matrix.n_cols; j++) {
            sub_matrix(i, j) = (*this)(i + 1, j + 1);
        }
    }
    return sub_matrix;
}

// namespace matrix

std::string matrix::read_file(std::string filename) {
    std::string file_string = " ";
    std::ifstream file;
    std::string new_line;
    file.open(filename);
    if (!file.is_open()) {
        throw std::invalid_argument("File does not exist.");
    }
    while (!file.eof()) {
        std::getline(file, new_line);
        for (size_t i = 0; i < new_line.length(); i++) {
            if (new_line[i] == ' ') {
                if (file_string.back() != ' ') {
                    file_string.push_back(' ');
                }
            } else {
                file_string.push_back(new_line[i]);
            }
        }
        if (file_string.back() != ' ') {
            file_string.push_back(' ');
        }
    }
    file.close();
    return file_string;
}

Matrix matrix::read_matrix_str(std::string matrix_string) {
    std::vector<int> string_values(0);
    std::string next_string = "";
    for (size_t i = 0; i < matrix_string.length(); i++) {
        if (matrix_string[i] == ' ') {
            if (next_string != "") {
                string_values.push_back(std::stoi(next_string));
                next_string = "";
            }
        } else {
            if (isdigit(matrix_string[i]) or matrix_string[i] == '-') {
                next_string.push_back(matrix_string[i]);
            } else {
                throw std::invalid_argument(
                    "String must have represent a sequence of integers.");
            }
        }
    }
    if (next_string != "") {
        string_values.push_back(std::stoi(next_string));
    }
    long unsigned int check_size = string_values[0] * string_values[1] + 2;
    if (string_values.size() != check_size) {
        throw std::invalid_argument(
            "Invalid string, first two integers misspecify the size.");
    }
    Matrix A(string_values[0], string_values[1]);
    for (int i = 0; i < string_values[0]; i++) {
        for (int j = 0; j < string_values[1]; j++) {
            A(i, j) = string_values[i * string_values[1] + j + 2];
        }
    }
    return A;
}

std::string matrix::write_matrix_str(Matrix A) {
    std::string output_string = "";
    for (int i = 0; i < A.n_rows; i++) {
        for (int j = 0; j < A.n_cols; j++) {
            std::string number_string = std::to_string(A(i, j));
            for (size_t k = 0; k < number_string.length(); k++) {
                output_string.push_back(number_string[k]);
            }
            if (j == A.n_cols - 1) {
                output_string.push_back('\n');
            } else {
                output_string.push_back(' ');
            }
        }
    }
    return output_string;
}

Matrix matrix::count_neighbours(Matrix &A) {
    /**
     * \brief Performs the necessary 3x3 2D convolution required for Game of
     * Life. \param A Input matrix.
     */

    int n_rows = A.n_rows;
    int n_cols = A.n_cols;

    // Row convolution
    Matrix B(n_rows, n_cols);
    for (int i = 0; i < n_rows; i++) {
        for (int j = 1; j < n_cols - 1; j++) {
            B(i, j) = A(i, j - 1) + A(i, j) + A(i, j + 1);
        }
    }

    // Column convolution
    Matrix C(n_rows, n_cols);
    for (int i = 1; i < n_rows - 1; i++) {
        for (int j = 1; j < n_cols - 1; j++) {
            C(i, j) = B(i - 1, j) + B(i, j) + B(i + 1, j);
        }
    }
    return C;
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
