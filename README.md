# Game of Life

## Description

Contains software used to simulate Conway's Game of Life using high performance computing techniques.

## How to use the project

Firstly, to compile the code run

`cmake -S . -B build`

`cmake --build build`

To compile the documentation run

`cmake --build build -t doxygen`

`make -C docs/latex`

The two main scripts used to demonstrate the code are `run_single` and `run_hybrid`, representing the single-thread solution and the threaded, domain decomposed solution respectively.

The former can be run with two or three passed arguments as follows:

`./bin/run_single <max_age> <seed.txt>`

`./bin/run_single <max_age> <n_rows> <n_cols>`

where the second format causes a random binary seed to be generated according to the input size. These option configurations are summarised as `<world_options>` in the documentation.

For the demonstration of the parallelised code in `run_hybrid`, the number of threads, MPI ranks, and topology can all be specified as follows:

`export OPM_NUM_THREADS=<...>`

`mpirun -n <n_ranks> ./bin/run_hybrid <n_rank_rows> <n_rank_cols> <world_options>`

The program should confirm the parallelisation parameters in its output. Note that the number of rows(columns) of cells does not have to be divisble by the number of rows(columns) of ranks, but that the total number of ranks in the specified topology must be equal to n_ranks.

For the unit testing suite use

`./bin/test_conway`

For the specifics of running the other scripts see the documentation or the comments at the top of each file.
