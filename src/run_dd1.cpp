/**
 * @brief run_dd1.cpp This file can be used to run the 1D domain decomposed
 * solution.
 * @details Can be run using `mpirun -n n_ranks ./bin/run_dd1 <world_options>`.
 */

#include <iostream>
#include <mpi.h>
#include <string>
#include <tuple>

#include "conway/include/matrix.h"
#include "conway/include/world.h"

int main(int argc, char **argv) {
    /**
     * Main procedure to run.
     */

    MPI_Init(&argc, &argv);
    int rank, nranks;
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    int max_age      = std::stoi(argv[1]);
    bool random_seed = (argc == 4);

    int n_rows_total;
    int n_cols_total;

    // Generate seed randomly or via input file.

    if (rank == 0) {
        if (random_seed) {
            n_rows_total = std::stoi(argv[2]);
            n_cols_total = std::stoi(argv[3]);
        } else {
            Matrix seed  = matrix::read_matrix_str(matrix::read_file(argv[2]));
            n_rows_total = seed.n_rows;
            n_cols_total = seed.n_cols;
        }
    }

    MPI_Bcast(&n_rows_total, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_cols_total, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Gen seed on rank 0
    auto full_data = new int[n_rows_total * n_cols_total];

    if (rank == 0) {
        Matrix seed = (random_seed)
                          ? matrix::generate_matrix(n_rows_total, n_cols_total)
                          : matrix::read_matrix_str(matrix::read_file(argv[2]));
        std::cout << "Initial seed:\n"
                  << matrix::write_matrix_str(seed) << std::endl;
        for (int i = 0; i < n_rows_total; i++) {
            for (int j = 0; j < n_cols_total; j++) {
                full_data[i * n_cols_total + j] = seed(i, j);
            }
        }
    }

    // Find chunk_sizes and offsets for Scatterv, Gatherv
    auto chunk_sizes = new int[nranks];
    auto offsets     = new int[nranks];

    for (int r = 0; r < nranks; r++) {
        int row_start, row_end;
        std::tie(row_start, row_end) =
            conway::divide_rows(n_rows_total, nranks, r);
        chunk_sizes[r] = (row_end - row_start) * n_cols_total;
        offsets[r]     = (r == 0) ? 0 : offsets[r - 1] + chunk_sizes[r - 1];
    }

    // Scatter seed
    int row_start, row_end;
    std::tie(row_start, row_end) =
        conway::divide_rows(n_rows_total, nranks, rank);
    auto chunk_data = new int[(row_end - row_start) * n_cols_total];

    MPI_Scatterv(full_data, chunk_sizes, offsets, MPI_INT, chunk_data,
                 (row_end - row_start) * n_cols_total, MPI_INT, 0,
                 MPI_COMM_WORLD);

    Matrix seed_chunk = Matrix(row_end - row_start, n_cols_total);
    for (int i = 0; i < (row_end - row_start); i++) {
        for (int j = 0; j < n_cols_total; j++) {
            seed_chunk(i, j) = chunk_data[i * (n_cols_total) + j];
        }
    }

    World WorldChunk(seed_chunk);

    for (int age = 0; age < max_age; age++) {
        // Update within rank periodicity first.
        WorldChunk.update_boundary();

        int rank_up   = (rank - 1 + nranks) % nranks;
        int rank_down = (rank + 1) % nranks;

        auto top_edge = new int[n_cols_total + 2];
        WorldChunk.read_edge_1d(top_edge, 0);

        auto bottom_edge = new int[n_cols_total + 2];
        WorldChunk.read_edge_1d(bottom_edge, 1);

        // Cycle edges
        MPI_Sendrecv_replace(top_edge, n_cols_total + 2, MPI_INT, rank_down, 42,
                             rank_up, 42, MPI_COMM_WORLD, &status);
        MPI_Sendrecv_replace(bottom_edge, n_cols_total + 2, MPI_INT, rank_up,
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
        for (int j = 0; j < n_cols_total; j++) {
            chunk_data[i * (n_cols_total) + j] = final_chunk(i, j);
        }
    }
    MPI_Gatherv(chunk_data, (row_end - row_start) * n_cols_total, MPI_INT,
                full_data, chunk_sizes, offsets, MPI_INT, 0, MPI_COMM_WORLD);

    // Write to a matrix and display via rank 0

    if (rank == 0) {
        Matrix final_state(n_rows_total, n_cols_total);
        for (int i = 0; i < n_rows_total; i++) {
            for (int j = 0; j < n_cols_total; j++) {
                final_state(i, j) = full_data[i * n_cols_total + j];
            }
        }
        std::cout << "Final state at age " << max_age << ":\n"
                  << matrix::write_matrix_str(final_state);
    }

    MPI_Finalize();
    return 0;
}
