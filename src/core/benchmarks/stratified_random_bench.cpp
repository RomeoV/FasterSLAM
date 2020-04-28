#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "stratified_random.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`


auto data_loader(size_t N, double* arr) {
    srand(0);
    for (int i = 0; i<N; i++) {
        arr[i] = 0.0;
    }
}

int main() {
    const size_t N = 100;
    double exact_di[N];

    double di[N];

    // Test:
    srand(0);
    stratified_random(N, di);

    srand(0); // reset seed!
    stratified_random_base(N, exact_di);

    for (int i = 0; i < N; i++) {
        double error = fabs(di[i]-exact_di[i]);
        expect(that % error < 1e-12) << i;
    }

    Benchmark<decltype(&stratified_random)> bench("pi_to_pi Benchmark");
    double work = 4*N+1; //How much does rand() cost??
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&stratified_random_base, "base", work);
    bench.add_function(&stratified_random, "active", work);

    bench.run_benchmark(N,di);

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