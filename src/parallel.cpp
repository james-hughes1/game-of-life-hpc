/**
 * \file main.cpp Main program file.
 */

#include <iostream>
#include <mpi.h>
#include <string>
#include <tuple>

#include "conway/include/matrix.h"
#include "conway/include/world.h"

std::tuple<int, int> divide_rows(int rows, int size, int rank) {
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

int main(int argc, char **argv) {
    /**
     * Main procedure to run.
     */

    MPI_Init(&argc, &argv);
    int rank, nranks;
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    int N_ROWS_TOTAL = 5;
    int N_COLS_TOTAL = 5;
    int MAX_AGE      = 1;

    // Gen seed on rank 0
    auto full_data = new int[N_ROWS_TOTAL * N_COLS_TOTAL];

    if (rank == 0) {
        Matrix seed = matrix::generate_matrix(N_ROWS_TOTAL, N_COLS_TOTAL);
        std::cout << "Initial seed:\n" << matrix::write_matrix_str(seed);
        for (int i = 0; i < N_ROWS_TOTAL; i++) {
            for (int j = 0; j < N_COLS_TOTAL; j++) {
                full_data[i * N_COLS_TOTAL + j] = seed(i, j);
            }
        }
    }

    // Find chunk_sizes and offsets for Scatterv, Gatherv
    auto chunk_sizes = new int[nranks];
    auto offsets     = new int[nranks];

    for (int r = 0; r < nranks; r++) {
        int row_start, row_end;
        std::tie(row_start, row_end) = divide_rows(N_ROWS_TOTAL, nranks, r);
        chunk_sizes[r]               = (row_end - row_start) * N_COLS_TOTAL;
        offsets[r] = (r == 0) ? 0 : offsets[r - 1] + chunk_sizes[r - 1];
    }

    // Scatter seed
    int row_start, row_end;
    std::tie(row_start, row_end) = divide_rows(N_ROWS_TOTAL, nranks, rank);
    auto chunk_data = new int[(row_end - row_start) * N_COLS_TOTAL];

    MPI_Scatterv(full_data, chunk_sizes, offsets, MPI_INT, chunk_data,
                 (row_end - row_start) * N_COLS_TOTAL, MPI_INT, 0,
                 MPI_COMM_WORLD);

    Matrix seed_chunk = Matrix(row_end - row_start, N_COLS_TOTAL);
    for (int i = 0; i < (row_end - row_start); i++) {
        for (int j = 0; j < N_COLS_TOTAL; j++) {
            seed_chunk(i, j) = chunk_data[i * (N_COLS_TOTAL) + j];
        }
    }

    World WorldChunk(seed_chunk);

    for (int age = 0; age < MAX_AGE; age++) {
        // Update within rank periodicity first.
        WorldChunk.update_boundary();

        int rank_up   = (rank - 1 + nranks) % nranks;
        int rank_down = (rank + 1) % nranks;

        auto top_edge = new int[N_COLS_TOTAL + 2];
        WorldChunk.read_edge(top_edge, 0);

        auto bottom_edge = new int[N_COLS_TOTAL + 2];
        WorldChunk.read_edge(bottom_edge, 1);

        // Cycle edges
        MPI_Sendrecv_replace(top_edge, N_COLS_TOTAL + 2, MPI_INT, rank_down, 42,
                             rank_up, 42, MPI_COMM_WORLD, &status);
        MPI_Sendrecv_replace(bottom_edge, N_COLS_TOTAL + 2, MPI_INT, rank_up,
                             42, rank_down, 42, MPI_COMM_WORLD, &status);

        // Write new boundaries from other ranks
        WorldChunk.write_edge(top_edge, 0);
        WorldChunk.write_edge(bottom_edge, 1);

        // Evaluate rules
        WorldChunk.evaluate_rules();
    }

    // Gather worlds to rank 0.
    Matrix final_chunk = WorldChunk.output_cells();
    for (int i = 0; i < (row_end - row_start); i++) {
        for (int j = 0; j < N_COLS_TOTAL; j++) {
            chunk_data[i * (N_COLS_TOTAL) + j] = final_chunk(i, j);
        }
    }
    MPI_Gatherv(chunk_data, (row_end - row_start) * N_COLS_TOTAL, MPI_INT,
                full_data, chunk_sizes, offsets, MPI_INT, 0, MPI_COMM_WORLD);

    // Write to a matrix and display via rank 0

    if (rank == 0) {
        Matrix final_state(N_ROWS_TOTAL, N_COLS_TOTAL);
        for (int i = 0; i < N_ROWS_TOTAL; i++) {
            for (int j = 0; j < N_COLS_TOTAL; j++) {
                final_state(i, j) = full_data[i * N_COLS_TOTAL + j];
            }
        }
        std::cout << "Final state at age " << MAX_AGE << ":\n"
                  << matrix::write_matrix_str(final_state);
    }

    MPI_Finalize();
    return 0;
}
