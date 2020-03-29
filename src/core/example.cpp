#include "example.hpp"

void add_two_arrays(double* lhs, double* rhs, double* res, size_t N) {
  for (size_t i = 0; i < N; i++) {
    res[i] = lhs[i] + rhs[i];
  }
}
