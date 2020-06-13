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

__m256d data_loader( __m256d a, __m256d b, __m256d c ) { 
    return a;
}

__m256d mmadd_2x2_reg( __m256d a, __m256d b, __m256d c ) {
    double A[4] __attribute__ ((aligned(32)));
    double B[4] __attribute__ ((aligned(32)));
    double C[4] __attribute__ ((aligned(32)));
    _mm256_store_pd(A, a);
    _mm256_store_pd(B, b);
    _mm256_store_pd(C, c);
    mmadd_2x2(A, B, C);
    return a; //dummy
}

int main() {

    // Test: 
    double A[4] __attribute__ ((aligned(32))) = {0., 1., 2., 3.};
    double B[4] __attribute__ ((aligned(32))) = {1., 2., 3., 4.};
    double C[4] __attribute__ ((aligned(32))) = {8., 9., 16., 21.};

    __m256d const a = _mm256_load_pd( A );
    __m256d const b = _mm256_load_pd( B );
    __m256d const c = _mm256_load_pd( C );

    Benchmark<decltype(&_mmadd_2x2_avx_v2)> bench("mmadd_2x2 benchmark");
    
    bench.data_loader = data_loader;
    
    bench.add_function(&mmadd_2x2_reg, "base", 0.0);
    bench.funcFlops[0] = mmadd_2x2_flops(A, B, C);
    bench.funcBytes[0] = 8*mmadd_2x2_memory(A, B, C);
    bench.add_function(&_mmadd_2x2_avx_v2, "avx_v2", 0.0);
    bench.funcFlops[1] = mmadd_2x2_flops(A, B, C);
    bench.funcBytes[1] = 8*mmadd_2x2_memory(A, B, C);

    bench.run_benchmark(a, b, c);

    return 0;
}

