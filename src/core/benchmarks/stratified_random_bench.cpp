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

    Benchmark<decltype(&stratified_random)> bench("strafified_random Benchmark");
    data_loader(N, di);
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&stratified_random_base, "base", 0.0);
    bench.funcFlops[0] = stratified_random_base_flops(N, di);
    bench.funcBytes[0] = 8*stratified_random_base_memory(N, di);
    bench.add_function(&stratified_random, "active", 0.0);
    bench.funcFlops[1] = stratified_random_flops(N, di);
    bench.funcBytes[1] = 8*stratified_random_memory(N, di);

    bench.run_benchmark(N, di);

    return 0;
}
