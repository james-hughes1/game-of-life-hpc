/**
 * \file main.cpp Main program file.
 */

#include <fstream>
#include <iostream>
#include <mpi.h>
#include <string>
#include <tuple>

#include "conway/include/matrix.h"
#include "conway/include/world.h"
#include "timing/include/timing.h"

int main(int argc, char **argv) {
    /**
     * Main procedure to run.
     */

    MPI_Init(&argc, &argv);
    int rank, nranks;
    MPI_Comm_size(MPI_COMM_WORLD, &nranks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    std::fstream file;
    double start_time = 0;

    if (rank == 0) {
        std::string output_filename =
            "prof/time_2d_decomp_" + std::to_string(nranks) + "rank.txt";
        file.open(output_filename);
    }

    for (int world_size = 100; world_size <= 4000; world_size += 100) {

        int n_rows_total = world_size;
        int n_cols_total = world_size;
        int MAX_AGE      = 50;

        // RANKS_ROWS * RANKS_COLS == nranks
        int RANKS_ROWS = 2;
        int RANKS_COLS = 2;

        // Split up rows and columns to enable proper ordering for full_data
        // into chunks.
        auto chunk_rows  = new int[RANKS_ROWS];
        auto offsets_row = new int[RANKS_ROWS];
        for (int i_chunk = 0; i_chunk < RANKS_ROWS; i_chunk++) {
            int row_start, row_end;
            std::tie(row_start, row_end) =
                conway::divide_rows(n_rows_total, RANKS_ROWS, i_chunk);
            chunk_rows[i_chunk]  = (row_end - row_start);
            offsets_row[i_chunk] = (i_chunk == 0) ? 0
                                                  : offsets_row[i_chunk - 1] +
                                                        chunk_rows[i_chunk - 1];
        }

        auto chunk_cols  = new int[RANKS_COLS];
        auto offsets_col = new int[RANKS_COLS];
        for (int j_chunk = 0; j_chunk < RANKS_COLS; j_chunk++) {
            int col_start, col_end;
            std::tie(col_start, col_end) =
                conway::divide_rows(n_cols_total, RANKS_COLS, j_chunk);
            chunk_cols[j_chunk]  = (col_end - col_start);
            offsets_col[j_chunk] = (j_chunk == 0) ? 0
                                                  : offsets_col[j_chunk - 1] +
                                                        chunk_cols[j_chunk - 1];
        }

        // Aggregate chunk dimensions into 'flat' sizes; same for offsets.
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
        auto full_data = new int[n_rows_total * n_cols_total];
        if (rank == 0) {
            Matrix seed = matrix::generate_matrix(n_rows_total, n_cols_total);
            int full_data_idx = 0;
            for (int i_chunk = 0; i_chunk < RANKS_ROWS; i_chunk++) {
                for (int j_chunk = 0; j_chunk < RANKS_COLS; j_chunk++) {
                    for (int i = 0; i < chunk_rows[i_chunk]; i++) {
                        for (int j = 0; j < chunk_cols[j_chunk]; j++) {
                            full_data[full_data_idx] =
                                seed(i + offsets_row[i_chunk],
                                     j + offsets_col[j_chunk]);
                            full_data_idx++;
                        }
                    }
                }
            }
        }

        // Create coords; note that these match rows and columns.
        int coord2d[2] = {rank / RANKS_COLS, rank % RANKS_COLS};

        // Assign markers to indicate which part of the world belongs to each
        // rank.
        int row_start = offsets_row[coord2d[0]];
        int col_start = offsets_col[coord2d[1]];
        int row_end   = ((coord2d[0] + 1) != RANKS_ROWS)
                            ? offsets_row[coord2d[0] + 1]
                            : n_rows_total;
        int col_end   = ((coord2d[1] + 1) != RANKS_COLS)
                            ? offsets_col[coord2d[1] + 1]
                            : n_cols_total;

        // Scatter world across ranks.
        auto chunk_data =
            new int[(row_end - row_start) * (col_end - col_start)];
        MPI_Scatterv(full_data, chunk_sizes, offsets, MPI_INT, chunk_data,
                     (row_end - row_start) * (col_end - col_start), MPI_INT, 0,
                     MPI_COMM_WORLD);

        // Abstract as a Matrix and then a World.
        Matrix seed_chunk = Matrix(row_end - row_start, col_end - col_start);
        for (int i = 0; i < (row_end - row_start); i++) {
            for (int j = 0; j < col_end - col_start; j++) {
                seed_chunk(i, j) = chunk_data[i * (col_end - col_start) + j];
            }
        }
        World WorldChunk(seed_chunk);

        // Find adjacent ranks
        int up = ((coord2d[0] + RANKS_ROWS - 1) % RANKS_ROWS) * RANKS_COLS +
                 coord2d[1];
        int down = ((coord2d[0] + 1) % RANKS_ROWS) * RANKS_COLS + coord2d[1];
        int left = coord2d[0] * RANKS_COLS +
                   ((coord2d[1] + RANKS_COLS - 1) % RANKS_COLS);
        int right = coord2d[0] * RANKS_COLS + ((coord2d[1] + 1) % RANKS_COLS);

        // Get diagonal neighbours
        int upright =
            ((coord2d[0] + RANKS_ROWS - 1) % RANKS_ROWS) * RANKS_COLS +
            ((coord2d[1] + 1) % RANKS_COLS);
        int downright = ((coord2d[0] + 1) % RANKS_ROWS) * RANKS_COLS +
                        ((coord2d[1] + 1) % RANKS_COLS);
        int downleft = ((coord2d[0] + 1) % RANKS_ROWS) * RANKS_COLS +
                       ((coord2d[1] + RANKS_COLS - 1) % RANKS_COLS);
        int upleft = ((coord2d[0] + RANKS_ROWS - 1) % RANKS_ROWS) * RANKS_COLS +
                     ((coord2d[1] + RANKS_COLS - 1) % RANKS_COLS);

        for (int age = 0; age < MAX_AGE; age++) {
            // Update within rank periodicity first.
            WorldChunk.update_boundary();

            // Cycle edges
            auto right_edge = new int[row_end - row_start];
            WorldChunk.read_edge_2d(right_edge, 3);
            MPI_Sendrecv_replace(right_edge, row_end - row_start, MPI_INT, left,
                                 42, right, 42, MPI_COMM_WORLD, &status);

            auto left_edge = new int[row_end - row_start];
            WorldChunk.read_edge_2d(left_edge, 1);
            MPI_Sendrecv_replace(left_edge, row_end - row_start, MPI_INT, right,
                                 42, left, 42, MPI_COMM_WORLD, &status);

            auto top_edge = new int[col_end - col_start];
            WorldChunk.read_edge_2d(top_edge, 0);
            MPI_Sendrecv_replace(top_edge, col_end - col_start, MPI_INT, down,
                                 42, up, 42, MPI_COMM_WORLD, &status);

            auto bottom_edge = new int[col_end - col_start];
            WorldChunk.read_edge_2d(bottom_edge, 2);
            MPI_Sendrecv_replace(bottom_edge, col_end - col_start, MPI_INT, up,
                                 42, down, 42, MPI_COMM_WORLD, &status);

            // Cycle vertices
            int top_left = WorldChunk.read_vertex_2d(0);
            MPI_Sendrecv_replace(&top_left, 1, MPI_INT, downright, 42, upleft,
                                 42, MPI_COMM_WORLD, &status);

            int bottom_left = WorldChunk.read_vertex_2d(1);
            MPI_Sendrecv_replace(&bottom_left, 1, MPI_INT, upright, 42,
                                 downleft, 42, MPI_COMM_WORLD, &status);

            int bottom_right = WorldChunk.read_vertex_2d(2);
            MPI_Sendrecv_replace(&bottom_right, 1, MPI_INT, upleft, 42,
                                 downright, 42, MPI_COMM_WORLD, &status);

            int top_right = WorldChunk.read_vertex_2d(3);
            MPI_Sendrecv_replace(&top_right, 1, MPI_INT, downleft, 42, upright,
                                 42, MPI_COMM_WORLD, &status);

            // Updates
            WorldChunk.write_edge_2d(right_edge, 3);
            WorldChunk.write_edge_2d(left_edge, 1);
            WorldChunk.write_edge_2d(top_edge, 0);
            WorldChunk.write_edge_2d(bottom_edge, 2);
            WorldChunk.write_vertex_2d(top_left, 0);
            WorldChunk.write_vertex_2d(bottom_left, 1);
            WorldChunk.write_vertex_2d(bottom_right, 2);
            WorldChunk.write_vertex_2d(top_right, 3);

            // Evaluate rules
            WorldChunk.evaluate_rules();
        }

        // Gather worlds to rank 0.
        Matrix final_chunk = WorldChunk.output_cells();
        for (int i = 0; i < (row_end - row_start); i++) {
            for (int j = 0; j < (col_end - col_start); j++) {
                chunk_data[i * (col_end - col_start) + j] = final_chunk(i, j);
            }
        }
        MPI_Gatherv(chunk_data, (row_end - row_start) * (col_end - col_start),
                    MPI_INT, full_data, chunk_sizes, offsets, MPI_INT, 0,
                    MPI_COMM_WORLD);

        // Read back from full_data
        if (rank == 0) {
            Matrix final_state(n_rows_total, n_cols_total);
            int full_data_idx = 0;
            for (int i_chunk = 0; i_chunk < RANKS_ROWS; i_chunk++) {
                for (int j_chunk = 0; j_chunk < RANKS_COLS; j_chunk++) {
                    for (int i = 0; i < chunk_rows[i_chunk]; i++) {
                        for (int j = 0; j < chunk_cols[j_chunk]; j++) {
                            final_state(i + offsets_row[i_chunk],
                                        j + offsets_col[j_chunk]) =
                                full_data[full_data_idx];
                            full_data_idx++;
                        }
                    }
                }
            }
            file << n_rows_total << " " << n_cols_total << " " << MAX_AGE << " "
                 << timing::get_split() - start_time << std::endl;
        }
    }

    if (rank == 0) {
        file.close();
    }

    MPI_Finalize();
    return 0;
}
