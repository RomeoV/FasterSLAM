#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "stratified_resample.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

void data_loader_bench_1(const double* w_, size_t N_w, double* Neff, size_t* keep) {
    for (int i = 0; i< N_w; i++) {
        *Neff = -1.0;
        keep[i] = N_w;
    }
}

int main() {

    const size_t N_w = 100;
    double w[N_w]; //unnormalized weights

    for (int i = 0; i< N_w; i++ ) {
        w[i] = 2*(i-1) +2;
    }
    double Neff; //Just written to this variable
    size_t keep[N_w]; //Just written to this variable
    data_loader_bench_1(w, N_w, &Neff, keep);
    
    size_t exact_keep[N_w]; //Just written to this variable
    double exact_Neff; //Just written to this variable
    data_loader_bench_1(w, N_w, &exact_Neff, exact_keep);

    // Test
    srand(0);
    stratified_resample(w, N_w, &Neff, keep);
    srand(0);
    stratified_resample_base(w, N_w, &exact_Neff, exact_keep);

    double error_Neff = fabs(Neff - exact_Neff);
    expect(that % error_Neff <= 1.0e-12);
    for (int i = 0; i < N_w; i++) {
        int error = keep[i]-exact_keep[i];
        expect(that % error == 0) << i;
    }

    double work = 4*N_w+1 + 5*N_w; //best

    Benchmark<decltype(&stratified_resample)> bench("stratified_resample Benchmark");
    
    bench.data_loader = data_loader_bench_1; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&stratified_resample_base, "base", work);
    bench.add_function(&stratified_resample, "active", work);

    // Linear weights
    bench.run_benchmark(w, N_w, &Neff, keep);

    //Uniform weights
    for (int i = 0; i< N_w; i++ ) {
        w[i] = 1.0/N_w;
    }

    bench.run_benchmark(w, N_w, &Neff, keep);

    //singular weights
    for (int i = 0; i< N_w; i++ ) {
        w[i] =0.0;
    }

    w[N_w-1] = 1.0;

    bench.run_benchmark(w, N_w, &Neff, keep);
    bench.destructor_output = false;
    bench.summary_long();
    /*
    // Alternative (much slower here, but nicer to look at. Generally useful if you want to average over a few inputs). Yields averages over all runs.
    
    Benchmark<decltype(&pi_to_pi)> multi_bench("pi_to_pi Benchmark");

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    multi_bench.add_function(&pi_to_pi, "pi_to_pi", 6);
    multi_bench.add_function(&pi_to_pi_fmod, "pi_to_pi_fmod", 6);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    for (int i = 0; i<N; i++) {
        // You could set the data_loader function here to generate new input. multi_bench.data_loader =&my_load_func_i...
        multi_bench.run_benchmark(angles[i]);
    }
    */

    return 0;
}