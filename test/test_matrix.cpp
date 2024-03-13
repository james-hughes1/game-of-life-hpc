#include <gtest/gtest.h>
#include <string>

#include "matrix.h"

std::string STR_FILE = " a b c d 1 2 3 4 X! Hello Wo_rld. ";
TEST(File, ReadFileValid) {
    std::string str_out = matrix::read_file("test/test_data/input_file_1.txt");
    EXPECT_EQ(STR_FILE, str_out)
        << "read_file doesn't parse valid files correctly.";
}

TEST(File, ReadFileInvalid) {
    EXPECT_THROW(matrix::read_file("test/test_data/input_file_2.txt"),
                 std::invalid_argument)
        << "Non-existent text file doesn't throw error.";
}

int populate_matrix(Matrix &matrix) {
    for (int i = 0; i < matrix.n_rows; i++) {
        for (int j = 0; j < matrix.n_cols; j++) {
            matrix(i, j) = matrix.n_cols * i + j + 1;
        }
    }
    return 0;
}

TEST(Matrix, MatrixEQValid) {
    Matrix MATRIX_1(3, 4);
    Matrix MATRIX_2(3, 4);
    populate_matrix(MATRIX_1);
    populate_matrix(MATRIX_2);
    EXPECT_TRUE(MATRIX_1 == MATRIX_2)
        << "Equal matrices are not found to be equal.";
}

TEST(Matrix, MatrixNEQValid) {
    Matrix MATRIX_1(3, 4);
    Matrix MATRIX_3(3, 4);
    populate_matrix(MATRIX_1);
    populate_matrix(MATRIX_3);
    MATRIX_3(2, 3) = -1000;
    EXPECT_TRUE(!(MATRIX_1 == MATRIX_3))
        << "Matrices of equal size but different entries are found to be "
           "equal.";
}

TEST(Matrix, MatrixNEQShapeValid) {
    Matrix MATRIX_1(3, 4);
    Matrix MATRIX_4(5, 5);
    populate_matrix(MATRIX_1);
    populate_matrix(MATRIX_4);
    EXPECT_TRUE(!(MATRIX_1 == MATRIX_4))
        << "Matrices of different sizes are found to be equal.";
}

TEST(Matrix, ReadMatrixValid) {
    Matrix MATRIX_1(3, 4);
    populate_matrix(MATRIX_1);
    MATRIX_1(2, 3)    = -1000;
    std::string STR_1 = "3 4 1 2 3 4 5 6 7 8 9 10 11 -1000";
    Matrix B          = matrix::read_matrix_str(STR_1);
    EXPECT_TRUE(matrix::read_matrix_str(STR_1) == MATRIX_1)
        << "Valid input string not correctly converted to matrix.";
}

TEST(Matrix, ReadMatrixInvalidShape) {
    // Shape wrong
    std::string STR_2 = "3 2 1 2 3 4 5 6 7 8 9 10 11 -1000";
    EXPECT_THROW(matrix::read_matrix_str(STR_2), std::invalid_argument)
        << "Invalid input string (bad shape) doesn't throw error.";
}

TEST(Matrix, ReadMatrixInvalidType) {
    // Contains anything other than spaces, digits and minus sign.
    std::string STR_3 = "4 3 1 2 3 4 5 6.1 7 8e-3 9 10 11 -1000";
    EXPECT_THROW(matrix::read_matrix_str(STR_3), std::invalid_argument)
        << "Invalid input string (bad chars) doesn't throw error.";
}

TEST(Matrix, WriteMatrixValid) {
    Matrix MATRIX_1(3, 4);
    populate_matrix(MATRIX_1);
    MATRIX_1(2, 3)    = -1000;
    std::string STR_1 = "1 2 3 4\n5 6 7 8\n9 10 11 -1000\n";
    EXPECT_EQ(matrix::write_matrix_str(MATRIX_1), STR_1)
        << "Valid input matrix not correctly converted to string.";
}

TEST(Matrix, GenerateBinaryMatrixValid) {
    Matrix A = matrix::generate_matrix(10, 20);
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 20; j++) {
            EXPECT_TRUE(A(i, j) == 0 or A(i, j) == 1)
                << "generate_matrix produced non-binary entries.";
        }
    }
}

TEST(Matrix, CountNeighboursValid) {
    Matrix MATRIX_5(4, 5);
    populate_matrix(MATRIX_5);
    Matrix Counts = matrix::count_neighbours(MATRIX_5);
    EXPECT_EQ(Counts.n_rows, 4);
    EXPECT_EQ(Counts.n_cols, 5);
    EXPECT_EQ(Counts(1, 1), 56)
        << "count_neighbours does not compute correct answers.";
    EXPECT_EQ(Counts(1, 2), 64)
        << "count_neighbours does not compute correct answers.";
    EXPECT_EQ(Counts(1, 3), 72)
        << "count_neighbours does not compute correct answers.";
    EXPECT_EQ(Counts(2, 1), 96)
        << "count_neighbours does not compute correct answers.";
    EXPECT_EQ(Counts(2, 2), 104)
        << "count_neighbours does not compute correct answers.";
    EXPECT_EQ(Counts(2, 3), 112)
        << "count_neighbours does not compute correct answers.";
}
