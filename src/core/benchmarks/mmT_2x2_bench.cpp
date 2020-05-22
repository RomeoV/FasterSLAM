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

auto data_loader(const double *A, const double *B, double *C) {
}

int main() {

    // Test: 
    double A[4] __attribute__ ((aligned(32))) = {0., 1., 2., 3.};
    double B[4] __attribute__ ((aligned(32))) = {1., 3., 2., 4.};
    double C[4] __attribute__ ((aligned(32))) = {0., 0., 0., 0.};
    double C1[4] __attribute__ ((aligned(32))) = {0., 0., 0., 0.};
    double C2[4] __attribute__ ((aligned(32))) = {0., 0., 0., 0.};
    double C3[4] __attribute__ ((aligned(32))) = {0., 0., 0., 0.};
    double D[4] __attribute__ ((aligned(32))) = {3., 4., 11., 16.};

    data_loader(A, B, C);
    mmT_2x2(A, B, C);
    
#ifdef __AVX2__
    data_loader(A, B, C1);
    mmT_2x2_avx_v1(A, B, C1);
    
    data_loader(A, B, C2);
    mmT_2x2_avx_v2(A, B, C2);
    
    data_loader(A, B, C3);
    mmT_2x2_avx_v3(A, B, C3);
#endif
    
    double error = 0.0;
    for (int i = 0; i < 4; i++) {
        error = fabs(  C[i] - D[i] ); expect(that % error < 1e-16) << i;
        error = fabs( C1[i] - D[i] ); expect(that % error < 1e-16) << i + 10;
        error = fabs( C2[i] - D[i] ); expect(that % error < 1e-16) << i + 20;
        error = fabs( C3[i] - D[i] ); expect(that % error < 1e-16) << i + 30;
    }
    
    Benchmark<decltype(&mmT_2x2)> bench("mmT_2x2 benchmark");
    data_loader(A, B, C);
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&mmT_2x2, "base", 0.0); 
    bench.funcFlops[0] = mmT_2x2_flops(A, B, C);
    bench.funcBytes[0] = 8*mmT_2x2_memory(A, B, C);
#ifdef __AVX2__
    bench.add_function(&mmT_2x2_avx_v1, "avx_v1", 0.0);
    bench.funcFlops[1] = mmT_2x2_flops(A, B, C);
    bench.funcBytes[1] = 8*mmT_2x2_memory(A, B, C);
    bench.add_function(&mmT_2x2_avx_v2, "avx_v2", 0.0);
    bench.funcFlops[2] = mmT_2x2_flops(A, B, C);
    bench.funcBytes[2] = 8*mmT_2x2_memory(A, B, C);
    bench.add_function(&mmT_2x2_avx_v3, "avx_v3", 0.0);
    bench.funcFlops[3] = mmT_2x2_flops(A, B, C);
    bench.funcBytes[3] = 8*mmT_2x2_memory(A, B, C);
#endif

    bench.run_benchmark(A, B, C);

    return 0;
}

