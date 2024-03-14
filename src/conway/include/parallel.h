#ifndef PARALLEL_H
#define PARALLEL_H

#include <tuple>

#include "matrix.h"
#include "world.h"

namespace parallel {
std::tuple<int, int> divide_rows(int rows, int size, int rank);
int decomposition_1d(int argc, char **argv);
}; // namespace parallel

#endif
