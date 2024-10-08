/**
 * @file test_dd1.cpp Tests that the 1D decomposition correctly evolves the
 * input_file_1.txt seed.
 * @details Can be run via `mpirun -n n_ranks ./bin/test_dd1`
 */

#include <iostream>
#include <mpi.h>
#include <string>
#include <tuple>

#include "matrix.h"
#include "world.h"

int main(int argc, char **argv) {
    /**
     * Main procedure to run.
     */

    MPI_Init(&argc, &argv);
    int rank, nranks;
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    int N_ROWS_TOTAL = 43;
    int N_COLS_TOTAL = 17;
    int MAX_AGE      = 50;

    // Gen seed on rank 0
    auto full_data = new int[N_ROWS_TOTAL * N_COLS_TOTAL];

    if (rank == 0) {
        std::string seed_str =
            matrix::read_file("test/test_data/input_file_2.txt");
        Matrix seed = matrix::read_matrix_str(seed_str);
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
        std::tie(row_start, row_end) =
            conway::divide_rows(N_ROWS_TOTAL, nranks, r);
        chunk_sizes[r] = (row_end - row_start) * N_COLS_TOTAL;
        offsets[r]     = (r == 0) ? 0 : offsets[r - 1] + chunk_sizes[r - 1];
    }

    // Scatter seed
    int row_start, row_end;
    std::tie(row_start, row_end) =
        conway::divide_rows(N_ROWS_TOTAL, nranks, rank);
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
        WorldChunk.read_edge_1d(top_edge, 0);

        auto bottom_edge = new int[N_COLS_TOTAL + 2];
        WorldChunk.read_edge_1d(bottom_edge, 1);

        // Cycle edges
        MPI_Sendrecv_replace(top_edge, N_COLS_TOTAL + 2, MPI_INT, rank_down, 42,
                             rank_up, 42, MPI_COMM_WORLD, &status);
        MPI_Sendrecv_replace(bottom_edge, N_COLS_TOTAL + 2, MPI_INT, rank_up,
                             42, rank_down, 42, MPI_COMM_WORLD, &status);

        // Write new boundaries from other ranks
        WorldChunk.write_edge_1d(top_edge, 0);
        WorldChunk.write_edge_1d(bottom_edge, 1);

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

    // Write to a matrix and test against expected.

    if (rank == 0) {
        Matrix final_state(N_ROWS_TOTAL, N_COLS_TOTAL);
        for (int i = 0; i < N_ROWS_TOTAL; i++) {
            for (int j = 0; j < N_COLS_TOTAL; j++) {
                final_state(i, j) = full_data[i * N_COLS_TOTAL + j];
            }
        }
        std::cout << "Final state at age " << MAX_AGE << ":\n"
                  << matrix::write_matrix_str(final_state);
        std::string expected_str =
            matrix::read_file("test/test_data/output_file_2.txt");
        Matrix expected_state = matrix::read_matrix_str(expected_str);
        if (final_state == expected_state) {
            std::cout << "Test passed." << std::endl;
        } else {
            std::cout << "Test failed." << std::endl;
        }
    }

    MPI_Finalize();
}
