/**
 * \file main.cpp Main program file.
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

    // MPI_Status status;

    int N_ROWS_TOTAL = 12;
    int N_COLS_TOTAL = 16;
    // int MAX_AGE      = 1;

    int RANKS_ROWS = 2;
    int RANKS_COLS = 2;

    // Split up rows and columns to enable proper ordering for full_data into
    // chunks.
    auto chunk_rows  = new int[RANKS_ROWS];
    auto offsets_row = new int[RANKS_ROWS];

    for (int i_chunk = 0; i_chunk < RANKS_ROWS; i_chunk++) {
        int row_start, row_end;
        std::tie(row_start, row_end) =
            conway::divide_rows(N_ROWS_TOTAL, RANKS_ROWS, i_chunk);
        chunk_rows[i_chunk] = (row_end - row_start);
        offsets_row[i_chunk] =
            (i_chunk == 0) ? 0
                           : offsets_row[i_chunk - 1] + chunk_rows[i_chunk - 1];
        if (rank == 0) {
            std::cout << chunk_rows[i_chunk] << " " << i_chunk << std::endl;
            std::cout << offsets_row[i_chunk] << " " << i_chunk << std::endl;
        }
    }

    auto chunk_cols  = new int[RANKS_COLS];
    auto offsets_col = new int[RANKS_COLS];

    for (int j_chunk = 0; j_chunk < RANKS_COLS; j_chunk++) {
        int col_start, col_end;
        std::tie(col_start, col_end) =
            conway::divide_rows(N_COLS_TOTAL, RANKS_COLS, j_chunk);
        chunk_cols[j_chunk] = (col_end - col_start);
        offsets_col[j_chunk] =
            (j_chunk == 0) ? 0
                           : offsets_col[j_chunk - 1] + chunk_cols[j_chunk - 1];
    }

    // Aggregate chunk dimensions into sizes; same for offsets.

    auto chunk_sizes = new int[nranks];
    auto offsets     = new int[nranks];

    for (int i_chunk = 0; i_chunk < RANKS_ROWS; i_chunk++) {
        for (int j_chunk = 0; j_chunk < RANKS_COLS; j_chunk++) {
            chunk_sizes[i_chunk * RANKS_COLS + j_chunk] =
                chunk_rows[i_chunk] * chunk_cols[j_chunk];
            offsets[i_chunk * RANKS_COLS + j_chunk] =
                (i_chunk + j_chunk == 0)
                    ? 0
                    : offsets[i_chunk * RANKS_COLS + j_chunk - 1] +
                          chunk_sizes[i_chunk * RANKS_COLS + j_chunk - 1];
        }
    }

    // Write send buffer for full_data from seed, on rank 0
    auto full_data = new int[N_ROWS_TOTAL * N_COLS_TOTAL];
    if (rank == 0) {
        Matrix seed = matrix::generate_matrix(N_ROWS_TOTAL, N_COLS_TOTAL);
        std::cout << "Initial seed:\n" << matrix::write_matrix_str(seed);
        int full_data_idx = 0;
        for (int i_chunk = 0; i_chunk < RANKS_ROWS; i_chunk++) {
            for (int j_chunk = 0; j_chunk < RANKS_COLS; j_chunk++) {
                for (int i = 0; i < chunk_rows[i_chunk]; i++) {
                    for (int j = 0; j < chunk_cols[j_chunk]; j++) {
                        full_data[full_data_idx] = seed(
                            i + offsets_row[i_chunk], j + offsets_col[j_chunk]);
                        full_data_idx++;
                    }
                }
            }
        }
    }

    // Scatter seed
    int row_start, row_end;
    std::tie(row_start, row_end) =
        conway::divide_rows(N_ROWS_TOTAL, RANKS_ROWS, rank);
    int col_start, col_end;
    std::tie(col_start, col_end) =
        conway::divide_rows(N_COLS_TOTAL, RANKS_COLS, rank);
    auto chunk_data = new int[(row_end - row_start) * (col_end - col_start)];

    MPI_Scatterv(full_data, chunk_sizes, offsets, MPI_INT, chunk_data,
                 (row_end - row_start) * (col_end - col_start), MPI_INT, 0,
                 MPI_COMM_WORLD);

    Matrix seed_chunk = Matrix(row_end - row_start, col_end - col_start);
    for (int i = 0; i < (row_end - row_start); i++) {
        for (int j = 0; j < col_end - col_start; j++) {
            seed_chunk(i, j) = chunk_data[i * (col_end - col_start) + j];
        }
    }

    std::cout << rank << std::endl
              << matrix::write_matrix_str(seed_chunk) << std::endl;

    //
    // World WorldChunk(seed_chunk);

    // for (int age = 0; age < MAX_AGE; age++) {
    //     // Update within rank periodicity first.
    //     WorldChunk.update_boundary();
    //
    //    int rank_up   = (rank - 1 + nranks) % nranks;
    //    int rank_down = (rank + 1) % nranks;
    //
    //    auto top_edge = new int[N_COLS_TOTAL + 2];
    //    WorldChunk.read_edge(top_edge, 0);
    //
    //    auto bottom_edge = new int[N_COLS_TOTAL + 2];
    //    WorldChunk.read_edge(bottom_edge, 1);
    //
    //    // Cycle edges
    //    MPI_Sendrecv_replace(top_edge, N_COLS_TOTAL + 2, MPI_INT, rank_down,
    //    42,
    //                         rank_up, 42, MPI_COMM_WORLD, &status);
    //    MPI_Sendrecv_replace(bottom_edge, N_COLS_TOTAL + 2, MPI_INT, rank_up,
    //                         42, rank_down, 42, MPI_COMM_WORLD, &status);
    //
    //    // Write new boundaries from other ranks
    //    WorldChunk.write_edge(top_edge, 0);
    //    WorldChunk.write_edge(bottom_edge, 1);
    //
    //    // Evaluate rules
    //    WorldChunk.evaluate_rules();
    //}
    //
    //// Gather worlds to rank 0.
    // Matrix final_chunk = WorldChunk.output_cells();
    // for (int i = 0; i < (row_end - row_start); i++) {
    //     for (int j = 0; j < N_COLS_TOTAL; j++) {
    //         chunk_data[i * (N_COLS_TOTAL) + j] = final_chunk(i, j);
    //     }
    // }
    // MPI_Gatherv(chunk_data, (row_end - row_start) * N_COLS_TOTAL, MPI_INT,
    //             full_data, chunk_sizes, offsets, MPI_INT, 0, MPI_COMM_WORLD);
    //
    //// Write to a matrix and display via rank 0
    //
    // if (rank == 0) {
    //    Matrix final_state(N_ROWS_TOTAL, N_COLS_TOTAL);
    //    for (int i = 0; i < N_ROWS_TOTAL; i++) {
    //        for (int j = 0; j < N_COLS_TOTAL; j++) {
    //            final_state(i, j) = full_data[i * N_COLS_TOTAL + j];
    //        }
    //    }
    //    std::cout << "Final state at age " << MAX_AGE << ":\n"
    //              << matrix::write_matrix_str(final_state);
    //}

    MPI_Finalize();
    return 0;
}
