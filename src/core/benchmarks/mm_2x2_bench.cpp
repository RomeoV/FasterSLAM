#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "linalg.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

__m256d data_loader( __m256d a, __m256d b ) {
    return a;
}

__m256d mm_2x2_reg( __m256d a, __m256d b ) {
    double A[4] __attribute__ ((aligned(32)));
    double B[4] __attribute__ ((aligned(32)));
    double C[4] __attribute__ ((aligned(32)));
    _mm256_store_pd(A, a);
    _mm256_store_pd(B, b);
    mm_2x2(A, B, C);
    return a; //dummy
}

int main() {

    // Test: 
    double A[4] __attribute__ ((aligned(32))) = {0., 1., 2., 3.};
    double B[4] __attribute__ ((aligned(32))) = {1., 2., 3., 4.};
    double C[4] __attribute__ ((aligned(32))) = {3., 4., 11., 16.};

    __m256d const a = _mm256_load_pd( A );
    __m256d const b = _mm256_load_pd( B );
 
    Benchmark<decltype(&_mm_2x2_avx_v1)> bench("mm_2x2 benchmark");

    bench.data_loader = data_loader; // To guarantee same inputs
    
    bench.add_function(&mm_2x2_reg, "base", 0.0);
    bench.funcFlops[0] = mm_2x2_flops(A, B, C);
    bench.funcBytes[0] = 8*mm_2x2_memory(A, B, C);
 
    bench.add_function(&_mm_2x2_avx_v1, "avx_v1", 0.0);
    bench.funcFlops[1] = mm_2x2_flops(A, B, C);
    bench.funcBytes[1] = 8*mm_2x2_memory(A, B, C);

    bench.run_benchmark(a, b);

    return 0;
}

