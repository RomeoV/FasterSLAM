#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"
#include "ut.hpp"

#if defined(__MACH__)
#include <stdlib.h>
#else 
#include <malloc.h>
#endif

#include <random>
#include <iostream>

using namespace boost::ut;

void create_array(int n, double& d) {
  double nums[n];
  d = nums[n-1];
}

void create_array_2(int n, double& d) {
  double *nums = (double *) malloc(n * sizeof(double));
  d = nums[n-1];
  free(nums);
}


int main(int argc, char* argv[]) {
  "VLA vs malloc"_test = [] {
    int n = rand()%1000 + 1;
    auto bench = ankerl::nanobench::Bench();
    bench.minEpochIterations(300);

    double sum = 0;
    bench.run("dynamic size", [&] {
        double d;
        create_array(n, d);
        sum += d;
        }
        ).doNotOptimizeAway(sum);

    bench.run("malloc", [&] {
        double d;
        create_array_2(n, d);
        sum += d;
        }
        ).doNotOptimizeAway(sum);
  };
}
