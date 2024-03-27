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

The two main scripts used to demonstrate the code are `main` and `hybrid`, representing the single-thread solution and the threaded, 2x2 domain decomposed solution respectively.

The former can be run with two or three passed arguments as follows:

`./bin/main <max_age> <seed.txt>`

`./bin/main <max_age> <n_rows> <n_cols>`

where the second format causes a random binary seed to be generated according to the inputted size.

For the demonstration of the parallelised code in `hybrid`, the number of threads, MPI ranks, and configuration can all be specified as follows:

`export OPM_NUM_THREADS=<...>`

`mpirun -n <n_ranks> ./bin/hybrid <n_rank_rows> <n_rank_cols> <max_age> <seed.txt>`, or

`mpirun -n <n_ranks> ./bin/hybrid <n_rank_rows> <n_rank_cols> <max_age> <n_rows> <n_cols>`

In either case the program should confirm the parallelisation parameters in its output.
Note that the number of rows(columns) of cells does not have to be divisble by the number of rows(columns) of ranks.

and for the test suite use

`./bin/test_conway`

For profiling (replace `<...>` as appropriate), run

`g++ -O3 -march=native -pg src/conway/matrix.cpp src/conway/world.cpp src/prof_<...>.cpp -o prof/prof_<...>`

`prof/prof_<...>`

`gprof prof/prof_<...> ./gmon.out > prof/prof_<...>.txt`

Similarly, for the timing script (and saving to .txt) run

`./bin/time_<...>`

or for the domain decompositions, run

`mpirun -n <1/2/4> ./bin/time_1d_decomp`



Either run as ./bin/executable max_age seed.txt
 or     ./bin/executable max_age n_rows n_cols

## Details
