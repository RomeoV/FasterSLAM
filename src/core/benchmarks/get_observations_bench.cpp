#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "get_observations.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

auto data_loader(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, int *idf, size_t *nidf, Vector2d z[]) {
    for (int i = 0; i < *nidf; i++) {
        idf[i] = i;
    }
}

int main() {
    // the function only calls base, so the bench is not necessary

    return 0;
}

