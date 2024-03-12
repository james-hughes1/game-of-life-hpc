#include <gtest/gtest.h>

#include "matrix.h"
#include "world.h"

TEST(World, Seed) {
    Matrix Cells_input(2, 1);
    Cells_input.zero();
    World NewWorld(Cells_input);
    EXPECT_TRUE(NewWorld.Cells_0.n_rows == 4 and NewWorld.Cells_0.n_cols == 3);
}

TEST(World, OutputCells) {
    Matrix Cells_input(2, 2);
    Cells_input(0, 0) = 1;
    Cells_input(0, 1) = 0;
    Cells_input(1, 0) = 0;
    Cells_input(1, 1) = 1;
    World NewWorld(Cells_input);

    Matrix Cells_output = NewWorld.output_cells();
    EXPECT_TRUE(Cells_output == Cells_input);
}

TEST(World, EvaluateRules) {
    Matrix Cells_input(4, 5);
    Cells_input.zero();
    Cells_input(0, 2) = 1;
    Cells_input(1, 2) = 1;
    Cells_input(1, 3) = 1;
    Cells_input(1, 4) = 1;
    Cells_input(2, 1) = 1;
    Cells_input(2, 4) = 1;

    Matrix Cells_output(4, 5);
    Matrix Cells_count(4, 5);
    Cells_count(1, 1) = 3;
    Cells_count(1, 2) = 3;
    Cells_count(1, 3) = 4;
    Cells_count(2, 1) = 1;
    Cells_count(2, 2) = 2;
    Cells_count(2, 3) = 4;

    conway::evaluate_rules(Cells_count, Cells_input, Cells_output);

    EXPECT_TRUE(Cells_count(1, 1) = 1);
    EXPECT_TRUE(Cells_count(1, 2) = 1);
    EXPECT_TRUE(Cells_count(1, 3) = 0);
    EXPECT_TRUE(Cells_count(2, 1) = 0);
    EXPECT_TRUE(Cells_count(2, 2) = 1);
    EXPECT_TRUE(Cells_count(2, 3) = 0);
}

TEST(World, UpdateBoundary) {
    Matrix Cells_input(4, 5);
    Cells_input.zero();
    Cells_input(0, 2) = 1;
    Cells_input(1, 2) = 1;
    Cells_input(1, 3) = 1;
    Cells_input(1, 4) = 1;
    Cells_input(2, 1) = 1;
    Cells_input(2, 4) = 1;

    conway::update_boundary(Cells_input);

    EXPECT_TRUE(Cells_input(0, 0) = 0);
    EXPECT_TRUE(Cells_input(0, 1) = 1);
    EXPECT_TRUE(Cells_input(0, 2) = 0);
    EXPECT_TRUE(Cells_input(0, 3) = 0);
    EXPECT_TRUE(Cells_input(0, 4) = 1);

    EXPECT_TRUE(Cells_input(1, 0) = 1);
    EXPECT_TRUE(Cells_input(2, 0) = 0);

    EXPECT_TRUE(Cells_input(1, 4) = 0);
    EXPECT_TRUE(Cells_input(2, 4) = 1);

    EXPECT_TRUE(Cells_input(3, 0) = 1);
    EXPECT_TRUE(Cells_input(3, 1) = 0);
    EXPECT_TRUE(Cells_input(3, 2) = 1);
    EXPECT_TRUE(Cells_input(3, 3) = 1);
    EXPECT_TRUE(Cells_input(3, 4) = 0);
}

TEST(World, Toad) {
    Matrix Cells_seed(100, 100);
    Cells_seed.zero();
    Cells_seed(49, 49) = 1;
    Cells_seed(49, 50) = 1;
    Cells_seed(49, 51) = 1;
    Cells_seed(50, 48) = 1;
    Cells_seed(50, 49) = 1;
    Cells_seed(50, 50) = 1;
    World ToadWorld(Cells_seed);
    for (int tick = 0; tick < 10; tick++) {
        ToadWorld.update_boundary();
        ToadWorld.evaluate_rules();
    }
    Matrix Cells_age_10 = ToadWorld.output_cells();
    EXPECT_TRUE(Cells_age_10 == Cells_seed);
    ToadWorld.update_boundary();
    ToadWorld.evaluate_rules();
    Matrix Cells_age_11 = ToadWorld.output_cells();
    EXPECT_TRUE(!(Cells_age_11 == Cells_seed));
}

TEST(World, Glider) { EXPECT_TRUE(true); }
