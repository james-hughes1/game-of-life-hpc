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
    Cells_count(2, 2) = 3;
    Cells_count(2, 3) = 4;

    conway::evaluate_rules(Cells_count, Cells_input, Cells_output);

    EXPECT_TRUE(Cells_output(1, 1) == 1);
    EXPECT_TRUE(Cells_output(1, 2) == 1);
    EXPECT_TRUE(Cells_output(1, 3) == 0);
    EXPECT_TRUE(Cells_output(2, 1) == 0);
    EXPECT_TRUE(Cells_output(2, 2) == 1);
    EXPECT_TRUE(Cells_output(2, 3) == 0);
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

    EXPECT_EQ(Cells_input(0, 0), 0);
    EXPECT_EQ(Cells_input(0, 1), 1);
    EXPECT_EQ(Cells_input(0, 2), 0);
    EXPECT_EQ(Cells_input(0, 3), 0);
    EXPECT_EQ(Cells_input(0, 4), 1);

    EXPECT_EQ(Cells_input(1, 0), 1);
    EXPECT_EQ(Cells_input(2, 0), 0);

    EXPECT_EQ(Cells_input(1, 4), 0);
    EXPECT_EQ(Cells_input(2, 4), 1);

    EXPECT_EQ(Cells_input(3, 0), 1);
    EXPECT_EQ(Cells_input(3, 1), 0);
    EXPECT_EQ(Cells_input(3, 2), 1);
    EXPECT_EQ(Cells_input(3, 3), 1);
    EXPECT_EQ(Cells_input(3, 4), 0);
}

TEST(World, Toad) {
    // https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life#References
    Matrix Cells_seed(20, 100);
    Cells_seed.zero();
    Cells_seed(9, 49)  = 1;
    Cells_seed(9, 50)  = 1;
    Cells_seed(9, 51)  = 1;
    Cells_seed(10, 48) = 1;
    Cells_seed(10, 49) = 1;
    Cells_seed(10, 50) = 1;
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

TEST(World, Glider) {
    // https://en.wikipedia.org/wiki/Conway%27s_Game_of_Life#References
    Matrix Cells_seed(20, 20);
    Cells_seed.zero();
    Cells_seed(8, 10)  = 1;
    Cells_seed(9, 10)  = 1;
    Cells_seed(10, 10) = 1;
    Cells_seed(10, 9)  = 1;
    Cells_seed(9, 8)   = 1;
    World GliderWorld(Cells_seed);
    for (int tick = 0; tick < 80; tick++) {
        // Every 4 ticks the glider moves one cell right, one down.
        // So it has periodicity 4*world_size if the world is square.
        GliderWorld.update_boundary();
        GliderWorld.evaluate_rules();
    }
    Matrix Cells_age_80 = GliderWorld.output_cells();
    EXPECT_TRUE(Cells_age_80 == Cells_seed);
    GliderWorld.update_boundary();
    GliderWorld.evaluate_rules();
    Matrix Cells_age_81 = GliderWorld.output_cells();
    EXPECT_TRUE(!(Cells_age_81 == Cells_seed));
}
