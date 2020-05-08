#include "fast_rand.h"
#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <math.h>

#define NR 10000

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

double sum = 0;

void sum_rands() {
    for (size_t i = 0; i < NR; i++) {
        sum += rand()/RAND_MAX;
    }
}

void sum_pseudo_rands() {
    for (size_t i = 0; i < NR; i++) {
        sum += xorshf96()/ulong_max;
    }
}

int main() {

    Benchmark<decltype(&sum_rands)> bench("RNG benchmark");
    double work = 1*NR;
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&sum_rands, "cmath", work);
    bench.add_function(&sum_pseudo_rands, "fast pseudo RNG", work);

    bench.run_benchmark();

    /*
    // Alternative (much slower here, but nicer to look at. Generally useful if you want to average over a few inputs). Yields averages over all runs.
    
    Benchmark<decltype(&pi_to_pi)> bench("pi_to_pi Benchmark");

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&pi_to_pi, "pi_to_pi", 6);
    bench.add_function(&pi_to_pi_fmod, "pi_to_pi_fmod", 6);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    for (int i = 0; i<N; i++) {
        // You could set the data_loader function here to generate new input. bench.data_loader =&my_load_func_i...
        bench.run_benchmark(angles[i]);
    }
    */

    return 0;
}

