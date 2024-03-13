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

To run the code use

`./bin/main`

and for the test suite use

`./bin/test_conway`

For profiling (replace `<...>` as appropriate), run

`g++ -pg src/conway/matrix.cpp src/conway/world.cpp src/prof_<...>.cpp -o prof/prof_<...>`

`prof/prof_<...>`

`gprof prof/prof_<...> ./gmon.out > prof/prof_<...>.txt`

Similarly, for the timing script (and saving to .txt) run

`./bin/time_simulation > prof/time_simulation.txt`

## Details
