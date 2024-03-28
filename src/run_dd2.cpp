/**
 * @brief run_dd2.cpp This file can be used to run the 2D domain decomposed
 * solution.
 * @details Can be run using `mpirun -n n_ranks ./bin/run_dd2 n_ranks_rows
 * n_cols_rows <world_options>`.
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

    // ranks_rows * ranks_cols == nranks
    // Defines the layout of the chunk topology.
    int ranks_rows = std::stoi(argv[1]);
    int ranks_cols = std::stoi(argv[2]);

    int max_age      = std::stoi(argv[3]);
    bool random_seed = (argc == 6);

    int n_rows_total;
    int n_cols_total;

    // Generate seed randomly or via input file.

    if (rank == 0) {
        if (random_seed) {
            n_rows_total = std::stoi(argv[4]);
            n_cols_total = std::stoi(argv[5]);
        } else {
            Matrix seed  = matrix::read_matrix_str(matrix::read_file(argv[4]));
            n_rows_total = seed.n_rows;
            n_cols_total = seed.n_cols;
        }
    }

    MPI_Bcast(&n_rows_total, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n_cols_total, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Split up rows and columns to enable proper ordering for full_data into
    // chunks.
    auto chunk_rows  = new int[ranks_rows];
    auto offsets_row = new int[ranks_rows];
    for (int i_chunk = 0; i_chunk < ranks_rows; i_chunk++) {
        int row_start, row_end;
        std::tie(row_start, row_end) =
            conway::divide_rows(n_rows_total, ranks_rows, i_chunk);
        chunk_rows[i_chunk] = (row_end - row_start);
        offsets_row[i_chunk] =
            (i_chunk == 0) ? 0
                           : offsets_row[i_chunk - 1] + chunk_rows[i_chunk - 1];
    }

    auto chunk_cols  = new int[ranks_cols];
    auto offsets_col = new int[ranks_cols];
    for (int j_chunk = 0; j_chunk < ranks_cols; j_chunk++) {
        int col_start, col_end;
        std::tie(col_start, col_end) =
            conway::divide_rows(n_cols_total, ranks_cols, j_chunk);
        chunk_cols[j_chunk] = (col_end - col_start);
        offsets_col[j_chunk] =
            (j_chunk == 0) ? 0
                           : offsets_col[j_chunk - 1] + chunk_cols[j_chunk - 1];
    }

    // Aggregate chunk dimensions into 'flat' sizes; same for offsets.
    auto chunk_sizes = new int[nranks];
    auto offsets     = new int[nranks];
    for (int i_chunk = 0; i_chunk < ranks_rows; i_chunk++) {
        for (int j_chunk = 0; j_chunk < ranks_cols; j_chunk++) {
            chunk_sizes[i_chunk * ranks_cols + j_chunk] =
                chunk_rows[i_chunk] * chunk_cols[j_chunk];
            offsets[i_chunk * ranks_cols + j_chunk] =
                (i_chunk + j_chunk == 0)
                    ? 0
                    : offsets[i_chunk * ranks_cols + j_chunk - 1] +
                          chunk_sizes[i_chunk * ranks_cols + j_chunk - 1];
        }
    }

    // Write send buffer for full_data from seed, on rank 0
    auto full_data = new int[n_rows_total * n_cols_total];
    if (rank == 0) {
        Matrix seed = (random_seed)
                          ? matrix::generate_matrix(n_rows_total, n_cols_total)
                          : matrix::read_matrix_str(matrix::read_file(argv[4]));
        std::cout << "Initial seed:\n"
                  << matrix::write_matrix_str(seed) << std::endl;
        int full_data_idx = 0;
        for (int i_chunk = 0; i_chunk < ranks_rows; i_chunk++) {
            for (int j_chunk = 0; j_chunk < ranks_cols; j_chunk++) {
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

    // Create coords; note that these match rows and columns.
    int coord2d[2] = {rank / ranks_cols, rank % ranks_cols};

    // Assign markers to indicate which part of the world belongs to each rank.
    int row_start = offsets_row[coord2d[0]];
    int col_start = offsets_col[coord2d[1]];
    int row_end = ((coord2d[0] + 1) != ranks_rows) ? offsets_row[coord2d[0] + 1]
                                                   : n_rows_total;
    int col_end = ((coord2d[1] + 1) != ranks_cols) ? offsets_col[coord2d[1] + 1]
                                                   : n_cols_total;

    // Scatter world across ranks.
    auto chunk_data = new int[(row_end - row_start) * (col_end - col_start)];
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
    int up =
        ((coord2d[0] + ranks_rows - 1) % ranks_rows) * ranks_cols + coord2d[1];
    int down = ((coord2d[0] + 1) % ranks_rows) * ranks_cols + coord2d[1];
    int left =
        coord2d[0] * ranks_cols + ((coord2d[1] + ranks_cols - 1) % ranks_cols);
    int right = coord2d[0] * ranks_cols + ((coord2d[1] + 1) % ranks_cols);

    // Get diagonal neighbours
    int upright = ((coord2d[0] + ranks_rows - 1) % ranks_rows) * ranks_cols +
                  ((coord2d[1] + 1) % ranks_cols);
    int downright = ((coord2d[0] + 1) % ranks_rows) * ranks_cols +
                    ((coord2d[1] + 1) % ranks_cols);
    int downleft = ((coord2d[0] + 1) % ranks_rows) * ranks_cols +
                   ((coord2d[1] + ranks_cols - 1) % ranks_cols);
    int upleft = ((coord2d[0] + ranks_rows - 1) % ranks_rows) * ranks_cols +
                 ((coord2d[1] + ranks_cols - 1) % ranks_cols);

    std::cout << "Rank " << rank << " topology coords: (" << coord2d[0] << ","
              << coord2d[1] << ")" << std::endl;

    for (int age = 0; age < max_age; age++) {
        // Update within rank periodicity first.
        WorldChunk.update_boundary();

        // Cycle edges
        auto right_edge = new int[row_end - row_start];
        WorldChunk.read_edge_2d(right_edge, 3);
        MPI_Sendrecv_replace(right_edge, row_end - row_start, MPI_INT, left, 42,
                             right, 42, MPI_COMM_WORLD, &status);

        auto left_edge = new int[row_end - row_start];
        WorldChunk.read_edge_2d(left_edge, 1);
        MPI_Sendrecv_replace(left_edge, row_end - row_start, MPI_INT, right, 42,
                             left, 42, MPI_COMM_WORLD, &status);

        auto top_edge = new int[col_end - col_start];
        WorldChunk.read_edge_2d(top_edge, 0);
        MPI_Sendrecv_replace(top_edge, col_end - col_start, MPI_INT, down, 42,
                             up, 42, MPI_COMM_WORLD, &status);

        auto bottom_edge = new int[col_end - col_start];
        WorldChunk.read_edge_2d(bottom_edge, 2);
        MPI_Sendrecv_replace(bottom_edge, col_end - col_start, MPI_INT, up, 42,
                             down, 42, MPI_COMM_WORLD, &status);

        // Cycle vertices
        int top_left = WorldChunk.read_vertex_2d(0);
        MPI_Sendrecv_replace(&top_left, 1, MPI_INT, downright, 42, upleft, 42,
                             MPI_COMM_WORLD, &status);

        int bottom_left = WorldChunk.read_vertex_2d(1);
        MPI_Sendrecv_replace(&bottom_left, 1, MPI_INT, upright, 42, downleft,
                             42, MPI_COMM_WORLD, &status);

        int bottom_right = WorldChunk.read_vertex_2d(2);
        MPI_Sendrecv_replace(&bottom_right, 1, MPI_INT, upleft, 42, downright,
                             42, MPI_COMM_WORLD, &status);

        int top_right = WorldChunk.read_vertex_2d(3);
        MPI_Sendrecv_replace(&top_right, 1, MPI_INT, downleft, 42, upright, 42,
                             MPI_COMM_WORLD, &status);

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
        for (int i_chunk = 0; i_chunk < ranks_rows; i_chunk++) {
            for (int j_chunk = 0; j_chunk < ranks_cols; j_chunk++) {
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
        std::cout << "Final state at age " << max_age << ":\n"
                  << matrix::write_matrix_str(final_state) << std::endl;
    }

    MPI_Finalize();
}
