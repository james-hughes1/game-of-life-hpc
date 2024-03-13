#ifndef MATRIX_H
#define MATRIX_H

#include <string>

/**
 * @brief Library for matrix functions
 * contained in the namespace "matrix"
 */

class Matrix {
  private:
    int *data;

  public:
    int n_rows;
    int n_cols;
    Matrix(int n_rows, int n_cols);
    int &operator()(int i, int j);
    ~Matrix();
    Matrix(const Matrix &old_matrix);
    Matrix &operator=(Matrix &old_matrix);
    bool operator==(Matrix &matrix_other);
    int zero();
    int write_sub_matrix(Matrix &sub_matrix);
    Matrix read_sub_matrix();
};

namespace matrix {
std::string read_file(std::string filename);
Matrix read_matrix_str(std::string matrix_string);
std::string write_matrix_str(Matrix A);
Matrix count_neighbours(Matrix A);
Matrix generate_matrix(int n_rows, int n_cols);
}; // namespace matrix

#endif
