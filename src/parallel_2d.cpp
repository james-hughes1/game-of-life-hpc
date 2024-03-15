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

    MPI_Status status;

    int N_ROWS_TOTAL = 4;
    int N_COLS_TOTAL = 6;
    int MAX_AGE      = 1;

    // RANKS_ROWS * RANKS_COLS == nranks
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
        std::cout << "Seed " << MAX_AGE << ":\n"
                  << matrix::write_matrix_str(seed);
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

    World WorldChunk(seed_chunk);

    MPI_Comm cartesian2d;
    int dims[2]    = {RANKS_ROWS, RANKS_COLS};
    int periods[2] = {1, 1};
    int reorder    = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cartesian2d);

    // Get the adjacent ranks
    // Note that North is right as looking at the true world (East is down).
    // To head north, increase second MPI topology coordinate.
    int west, east, north, south;
    MPI_Cart_shift(cartesian2d, 0, 1, &west, &east);
    MPI_Cart_shift(cartesian2d, 1, 1, &north, &south);

    // Get diagonal neighbours
    int coord2d[2];
    MPI_Cart_coords(cartesian2d, rank, 2, coord2d);
    const int coord_ne[2] = {coord2d[0] + 1, coord2d[1] + 1};
    const int coord_se[2] = {coord2d[0] + 1, coord2d[1] - 1};
    const int coord_sw[2] = {coord2d[0] - 1, coord2d[1] - 1};
    const int coord_nw[2] = {coord2d[0] - 1, coord2d[1] + 1};
    int northeast, southeast, southwest, northwest;
    MPI_Cart_rank(cartesian2d, coord_ne, &northeast);
    MPI_Cart_rank(cartesian2d, coord_se, &southeast);
    MPI_Cart_rank(cartesian2d, coord_sw, &southwest);
    MPI_Cart_rank(cartesian2d, coord_nw, &northwest);

    for (int age = 0; age < MAX_AGE; age++) {
        // Update within rank periodicity first.
        WorldChunk.update_boundary();

        std::cout << rank << std::endl
                  << matrix::write_matrix_str(WorldChunk.Cells_0) << std::endl;

        // Cycle edges
        auto right_edge = new int[row_end - row_start];
        WorldChunk.read_edge_2d(right_edge, 3);
        MPI_Sendrecv_replace(right_edge, row_end - row_start, MPI_INT, south,
                             42, north, 42, MPI_COMM_WORLD, &status);

        auto left_edge = new int[row_end - row_start];
        WorldChunk.read_edge_2d(left_edge, 1);
        MPI_Sendrecv_replace(left_edge, row_end - row_start, MPI_INT, north, 42,
                             south, 42, MPI_COMM_WORLD, &status);

        auto top_edge = new int[col_end - col_start];
        WorldChunk.read_edge_2d(top_edge, 0);
        MPI_Sendrecv_replace(top_edge, col_end - col_start, MPI_INT, east, 42,
                             west, 42, MPI_COMM_WORLD, &status);

        auto bottom_edge = new int[col_end - col_start];
        WorldChunk.read_edge_2d(bottom_edge, 2);
        MPI_Sendrecv_replace(bottom_edge, col_end - col_start, MPI_INT, west,
                             42, east, 42, MPI_COMM_WORLD, &status);

        // Cycle vertices

        int top_left = WorldChunk.read_vertex_2d(0);
        MPI_Sendrecv_replace(&top_left, 1, MPI_INT, northeast, 42, southwest,
                             42, MPI_COMM_WORLD, &status);

        int bottom_left = WorldChunk.read_vertex_2d(1);
        MPI_Sendrecv_replace(&bottom_left, 1, MPI_INT, northwest, 42, southeast,
                             42, MPI_COMM_WORLD, &status);

        int bottom_right = WorldChunk.read_vertex_2d(2);
        std::cout << rank << " : " << bottom_right << std::endl;
        std::cout << rank << " : " << southwest << std::endl;
        MPI_Sendrecv_replace(&bottom_right, 1, MPI_INT, southwest, 42,
                             northeast, 42, MPI_COMM_WORLD, &status);
        std::cout << rank << " : " << bottom_right << std::endl;

        int top_right = WorldChunk.read_vertex_2d(3);
        MPI_Sendrecv_replace(&top_right, 1, MPI_INT, southeast, 42, northwest,
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

        std::cout << rank << std::endl
                  << matrix::write_matrix_str(WorldChunk.Cells_0) << std::endl;

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
        Matrix final_state(N_ROWS_TOTAL, N_COLS_TOTAL);
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

        std::cout << "Final state at age " << MAX_AGE << ":\n"
                  << matrix::write_matrix_str(final_state);
    }

    MPI_Finalize();
    return 0;
}
