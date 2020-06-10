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

//
// DataLoader for Benchmark 1
//
auto data_loader(__m256d const r0,
                 __m256d const r1,
                 __m256d const r2,
                 __m256d const r3,
                 __m256d *inv0,
                 __m256d *inv1,
                 __m256d *inv2,
                 __m256d *inv3) {

}


void inv_2x2_x4(__m256d const r0,
                __m256d const r1,
                __m256d const r2,
                __m256d const r3,
                __m256d *inv0,
                __m256d *inv1,
                __m256d *inv2,
                __m256d *inv3) {

    double A[16], Ainv[16];
    _mm256_store_pd(A+0,  r0);
    _mm256_store_pd(A+4,  r1);
    _mm256_store_pd(A+8,  r2);
    _mm256_store_pd(A+12, r3);
    inv_2x2(A+0,  Ainv+0); 
    inv_2x2(A+4,  Ainv+4); 
    inv_2x2(A+8,  Ainv+8); 
    inv_2x2(A+12, Ainv+12); 
}

int main() {

    // Declaration of inputs + test
    double A[16] __attribute__ ((aligned(32)));
    for (int i = 0; i < 16; i++) {
        A[i] = i;
    }
    __m256d const r0 = _mm256_load_pd( A+0 );
    __m256d const r1 = _mm256_load_pd( A+4 );
    __m256d const r2 = _mm256_load_pd( A+8 );
    __m256d const r3 = _mm256_load_pd( A+12 );

    double Inv[16] __attribute__ ((aligned(32)));
    __m256d inv0, inv1, inv2, inv3;
    batch_inverse_2x2(r0, r1, r2, r3, &inv0, &inv1, &inv2, &inv3);

    _mm256_store_pd( Inv+0, inv0 );
    _mm256_store_pd( Inv+4, inv1 );
    _mm256_store_pd( Inv+8, inv2 );
    _mm256_store_pd( Inv+12, inv3 );

    double AInv[16]; 
    for (size_t i = 0; i < 4; i++) {
        inv_2x2(A+i*4, AInv+i*4);
        for (size_t j = 0; j < 4; j++) {
            expect(fabs(Inv[i*4+j] - AInv[i*4+j]) < 1e-10) << Inv[i*4+j] << " != " << AInv[i*4+j]; 
        }
    } 

    // ----------------------------------------- //
    // --------------- BENCHMARK --------------- //
    // ----------------------------------------- //

    Benchmark<decltype(&batch_inverse_2x2)> bench_u4avx("batch_inverse_2x2 benchmark (UNROLLEDx4)");
    bench_u4avx.data_loader = data_loader;

    bench_u4avx.add_function(&inv_2x2_x4, "base_x4", 0.0);
    bench_u4avx.funcFlops[0] = 0;
    bench_u4avx.funcBytes[0] = 0;
    
    bench_u4avx.add_function(&batch_inverse_2x2, "batch_inverse_2x2", 0.0);
    bench_u4avx.funcFlops[1] = 0;
    bench_u4avx.funcBytes[1] = 0;

    bench_u4avx.run_benchmark(r0, r1, r2, r3, &inv0, &inv1, &inv2, &inv3); 

    return 0;
}

