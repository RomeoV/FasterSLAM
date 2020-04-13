#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"

#include <malloc.h>
#include <iostream>

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
  int n;
  if (argc == 1)
    std::cin >> n;
  else
    n = std::stoi(argv[1]);
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
}
