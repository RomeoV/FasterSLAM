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
    double B[4] __attribute__ ((aligned(32))) = {1., 2., 3., 4.};
    double C[4] __attribute__ ((aligned(32))) = {0., 0., 0., 0.};
    double C1[4] __attribute__ ((aligned(32))) = {0., 0., 0., 0.};
    double C2[4] __attribute__ ((aligned(32))) = {0., 0., 0., 0.};
    double C3[4] __attribute__ ((aligned(32))) = {0., 0., 0., 0.};
    double D[4] __attribute__ ((aligned(32))) = {3., 4., 11., 16.};

    data_loader(A, B, C);
    mm_2x2(A, B, C);
    
#ifdef __AVX2__
    data_loader(A, B, C1);
    mm_2x2_avx_v1(A, B, C1);
    
    data_loader(A, B, C2);
    mm_2x2_avx_v2(A, B, C2);
    
    data_loader(A, B, C3);
    mm_2x2_avx_v3(A, B, C3);
#endif
    
    double error = 0.0;
    for (int i = 0; i < 4; i++) {
        error = fabs(  C[i] - D[i] ); expect(that % error < 1e-16) << i;
        error = fabs( C1[i] - D[i] ); expect(that % error < 1e-16) << i + 10;
        error = fabs( C2[i] - D[i] ); expect(that % error < 1e-16) << i + 20;
        error = fabs( C3[i] - D[i] ); expect(that % error < 1e-16) << i + 30;
    }
    
    Benchmark<decltype(&mm_2x2)> bench("mm_2x2 benchmark");
    double work = 12; 
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&mm_2x2, "base", work);
#ifdef __AVX2__
    bench.add_function(&mm_2x2_avx_v1, "avx_v1", work);
    bench.add_function(&mm_2x2_avx_v2, "avx_v2", work);
    bench.add_function(&mm_2x2_avx_v3, "avx_v2", work);
#endif

    bench.run_benchmark(A, B, C);

    return 0;
}

