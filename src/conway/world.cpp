/**
 * @file world.cpp Defines the World class for simulating Game of Life,
 * as well as some useful routines for the updating and communicating the halos
 * during domain decomposition.
 */

#include <iostream>
#include <random>
#include <string>

#include "include/matrix.h"
#include "include/world.h"

World::World(Matrix seed)
    /**
     * @brief Constructor for the World class
     * @details Creates two matrix members of the size of the seed parameter
     * plus a halo each. Fills the inner entries of one of these members with
     * the seed entries. Stores an age attribute to keep track of the progress
     * of the simulation, to keep track of which of the two member matrices
     * holds the current cell values.
     * @param seed The seed of values used to initialise the cell states.
     */
    : Cells_0(seed.n_rows + 2, seed.n_cols + 2),
      Cells_1(seed.n_rows + 2, seed.n_cols + 2) {
    Cells_0.write_sub_matrix(seed);
    n_rows = seed.n_rows;
    n_cols = seed.n_cols;
    age    = 0;
}

int World::evaluate_rules() {
    /**
     * @brief Updates the current cell values, aging the simulation by 1 tick.
     * @details Assumes that the ghost cells at the edge of the cells matrix has
     * been updated according to the desired boundary condition.
     */
    if (age % 2 == 0) {
        Matrix Cells_count = matrix::count_neighbours(Cells_0);
        conway::evaluate_rules(Cells_count, Cells_0, Cells_1);
    } else {
        Matrix Cells_count = matrix::count_neighbours(Cells_1);
        conway::evaluate_rules(Cells_count, Cells_1, Cells_0);
    }
    age += 1;
    return 0;
}

int World::update_boundary() {
    /**
     * @brief Invokes the matrix::update_boundary routine on the correct cell
     * matrix according to the parity of age.
     */
    if (age % 2 == 0) {
        conway::update_boundary(Cells_0);
    } else {
        conway::update_boundary(Cells_1);
    }
    return 0;
}

Matrix World::output_cells() {
    /**
     * @brief Converts the current cell states to a string.
     * @returns std::string the string of current cell states, formatted.
     */
    if (age % 2 == 0) {
        return Cells_0.read_sub_matrix();
    } else {
        return Cells_1.read_sub_matrix();
    }
}

int World::display_world() {
    /**
     * @brief Reports key information on the current simulation state, including
     * the current age and cell states.
     */
    std::cout << "World has age: " << age << " ticks." << std::endl;
    std::cout << "------------------------------" << std::endl;
    std::string output = matrix::write_matrix_str((*this).output_cells());
    std::cout << output << std::endl;
    return 0;
}

int World::write_edge_1d(int *edge, int loc) {
    /**
     * @brief Used to update edges (including vertices) during 1D domain
     * decomposition.
     * @param edge Edge buffer
     * @param loc Integer specifying the location of the edge, with 0 meaning
     * top and 1 meaning bottom.
     */
    if (age % 2 == 0) {
        for (int j = 0; j < Cells_0.n_cols; j++) {
            Cells_0(loc * (Cells_0.n_rows - 1), j) = edge[j];
        }
    } else {
        for (int j = 0; j < Cells_1.n_cols; j++) {
            Cells_1(loc * (Cells_1.n_rows - 1), j) = edge[j];
        }
    }
    return 0;
}

int World::read_edge_1d(int *edge, int loc) {
    /**
     * @brief Used to read edges (including vertices) during 1D domain
     * decomposition.
     * @param edge Edge buffer, which is updated in-place
     * @param loc Integer specifying the location of the edge, with 0 meaning
     * top and 1 meaning bottom.
     */
    if (age % 2 == 0) {
        for (int j = 0; j < Cells_0.n_cols; j++) {
            edge[j] = Cells_0(loc * (Cells_0.n_rows - 1), j);
        }
    } else {
        for (int j = 0; j < Cells_1.n_cols; j++) {
            edge[j] = Cells_1(loc * (Cells_1.n_rows - 1), j);
        }
    }
    return 0;
}

int World::write_edge_2d(int *edge, int loc) {
    /**
     * @brief Used to update edges (excluding vertices) during 2D domain
     * decomposition.
     * @param edge Edge buffer
     * @param loc Integer specifying the location of the edge, where 0=top,
     * 1=left, 2=bottom, 3=right
     */
    if (age % 2 == 0) {
        if (loc % 2 == 0) {
            for (int j = 1; j < n_cols + 1; j++) {
                Cells_0((loc / 2) * (n_rows + 1), j) = edge[j - 1];
            }
        } else {
            for (int i = 1; i < n_rows + 1; i++) {
                Cells_0(i, (loc / 2) * (n_cols + 1)) = edge[i - 1];
            }
        }
    } else {
        if (loc % 2 == 0) {
            for (int j = 1; j < n_cols + 1; j++) {
                Cells_1((loc / 2) * (n_rows + 1), j) = edge[j - 1];
            }
        } else {
            for (int i = 1; i < n_rows + 1; i++) {
                Cells_1(i, (loc / 2) * (n_cols + 1)) = edge[i - 1];
            }
        }
    }
    return 0;
}

int World::read_edge_2d(int *edge, int loc) {
    /**
     * @brief Used to read edges (excluding vertices) during 2D domain
     * decomposition.
     * @param edge Edge buffer, updated in-place
     * @param loc Integer specifying the location of the edge, where 0=top,
     * 1=left, 2=bottom, 3=right
     */
    if (age % 2 == 0) {
        if (loc % 2 == 0) {
            for (int j = 1; j < n_cols + 1; j++) {
                edge[j - 1] = Cells_0((loc / 2) * (n_rows + 1), j);
            }
        } else {
            for (int i = 1; i < n_rows + 1; i++) {
                edge[i - 1] = Cells_0(i, (loc / 2) * (n_cols + 1));
            }
        }
    } else {
        if (loc % 2 == 0) {
            for (int j = 1; j < n_cols + 1; j++) {
                edge[j - 1] = Cells_1((loc / 2) * (n_rows + 1), j);
            }
        } else {
            for (int i = 1; i < n_rows + 1; i++) {
                edge[i - 1] = Cells_1(i, (loc / 2) * (n_cols + 1));
            }
        }
    }
    return 0;
}
int World::read_vertex_2d(int loc) {
    /**
     * @brief Used to update vertices during 2D domain decomposition.
     * @param loc Integer specifying the location of the vertex, where
     * 0=top-left, 1=bottom-left, 2=bottom-right, 3=top-right
     * @returns vertex Vertex value
     */
    int vertex;
    if (loc == 0) {
        vertex = (age % 2 == 0) ? Cells_0(0, 0) : Cells_1(0, 0);
    } else if (loc == 1) {
        vertex =
            (age % 2 == 0) ? Cells_0(n_rows + 1, 0) : Cells_1(n_rows + 1, 0);
    } else if (loc == 2) {
        vertex = (age % 2 == 0) ? Cells_0(n_rows + 1, n_cols + 1)
                                : Cells_1(n_rows + 1, n_cols + 1);
    } else {
        vertex =
            (age % 2 == 0) ? Cells_0(0, n_cols + 1) : Cells_1(0, n_cols + 1);
    }
    return vertex;
}
int World::write_vertex_2d(int vertex, int loc) {
    /**
     * @brief Used to update vertices during 2D domain decomposition.
     * @param vertex Vertex value to be written
     * @param loc Integer specifying the location of the vertex, where
     * 0=top-left, 1=bottom-left, 2=bottom-right, 3=top-right
     */
    if (age % 2 == 0) {
        if (loc == 0) {
            Cells_0(0, 0) = vertex;
        } else if (loc == 1) {
            Cells_0(n_rows + 1, 0) = vertex;
        } else if (loc == 2) {
            Cells_0(n_rows + 1, n_cols + 1) = vertex;
        } else if (loc == 3) {
            Cells_0(0, n_cols + 1) = vertex;
        }
    } else {
        if (loc == 0) {
            Cells_1(0, 0) = vertex;
        } else if (loc == 1) {
            Cells_1(n_rows + 1, 0) = vertex;
        } else if (loc == 2) {
            Cells_1(n_rows + 1, n_cols + 1) = vertex;
        } else if (loc == 3) {
            Cells_1(0, n_cols + 1) = vertex;
        }
    }
    return 0;
}

int conway::evaluate_rules(Matrix &Cells_count, Matrix &Cells_current,
                           Matrix &Cells_next) {
    /**
     * @brief Evaluates the update of cell states according to Game of Life
     * rules.
     * @param Cells_count Matrix giving the relevant neighbouring cell counts
     * for each entry, including each cell's own value
     * @param Cells_current Matrix giving the current cell states
     * @param Cells_next Matrix giving the updated cell states, altered in-place
     */
    for (int i = 0; i < Cells_count.n_rows; i++) {
        for (int j = 0; j < Cells_count.n_cols; j++) {
            Cells_next(i, j) =
                (Cells_count(i, j) == 3) ||
                ((Cells_count(i, j) == 4) && (Cells_current(i, j) == 1));
        }
    }
    return 0;
}

int conway::update_boundary(Matrix &Cells_current) {
    /**
     * @brief Updates the ghost cells on the perimeter according to the periodic
     * boundary condition.
     * @param Cells_current Matrix to be updated, in-place.
     */
    // Vertices
    Cells_current(Cells_current.n_rows - 1, Cells_current.n_cols - 1) =
        Cells_current(1, 1);
    Cells_current(0, Cells_current.n_cols - 1) =
        Cells_current(Cells_current.n_rows - 2, 1);
    Cells_current(0, 0) =
        Cells_current(Cells_current.n_rows - 2, Cells_current.n_cols - 2);
    Cells_current(Cells_current.n_rows - 1, 0) =
        Cells_current(1, Cells_current.n_cols - 2);
    // Edges
    for (int j = 1; j < Cells_current.n_cols - 1; j++) {
        Cells_current(0, j) = Cells_current(Cells_current.n_rows - 2, j);
        Cells_current(Cells_current.n_rows - 1, j) = Cells_current(1, j);
    }
    for (int i = 1; i < Cells_current.n_rows - 1; i++) {
        Cells_current(i, 0) = Cells_current(i, Cells_current.n_cols - 2);
        Cells_current(i, Cells_current.n_cols - 1) = Cells_current(i, 1);
    }
    return 0;
}

std::tuple<int, int> conway::divide_rows(int rows, int size, int rank) {
    /**
     * @brief Routine used to divide rows (or columns) across a number of ranks
     * evenly that isn't necessarily a divisor.
     * @details This has been adapted from the `divide_rows` routine found in
     * https://gitlab.developers.cam.ac.uk/phy/data-intensive-science-mphil/
     * c2_advanced_research_computing/-/tree/main/Snippets/MPI?ref_type=heads
     * @param rows The total number of rows to be distributed
     * @param size The number of ranks
     * @param rank The index of the specific rank to be allocated
     * @returns std::make_tuple Tuple of indices of the first and last row to be
     * allocated to the rank.
     */
    int row_rank  = rows / size;
    int auxrow    = rows % size;
    int start_row = rank * row_rank;
    int end_row   = (rank + 1) * row_rank;

    // Allocate remainder across rows
    if (auxrow != 0) {
        if (rank < auxrow) {
            start_row = start_row + rank;
            end_row   = end_row + rank + 1;
        } else {
            start_row = start_row + auxrow;
            end_row   = end_row + auxrow;
        }
    }
    return std::make_tuple(start_row, end_row);
}
