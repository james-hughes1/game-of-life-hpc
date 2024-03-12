#include "matrix.h"
#include "world.h"
#include <gtest/gtest.h>

TEST(Matrix, GenerateBinaryMatrix) {
    Matrix A = matrix::generate_matrix(3, 4);
    EXPECT_TRUE(A(2, 3) == 0 or A(2, 3) == 1)
        << "generate_matrix produced non-binary entries.";
}
