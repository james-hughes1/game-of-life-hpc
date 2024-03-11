#ifndef MATRIX_H
#define MATRIX_H

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
};

namespace matrix {
int display_matrix(Matrix A);
Matrix count_neighbours(Matrix A);
Matrix generate_matrix(int n_rows, int n_cols);
}; // namespace matrix

#endif
