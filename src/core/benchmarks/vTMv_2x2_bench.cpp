#include "linalg.h"
#include "rdtsc_benchmark.h"
#include "typedefs.h"
#include <stdlib.h>

#define NR 300

void fill(double* Ms, double* Vs, double* results)
{
    fill_rand(Ms, 4*NR, -10, 10);
    fill_rand(Vs, 2*NR, -10, 10);
}

void vTMv_base(double* Ms, double* Vs, double* results) {
    for (size_t i = 0; i < NR; i++) {
        Matrix2d M_vT;
        mul(Ms + (4*i), Vs + (2*i), 2, 2, 1, M_vT);
        mul(Vs + (2*i), M_vT, 1, 2, 1, results + i);
    }
}

#ifdef __AVX2__
void vTMv_avx2(double* Ms, double* Vs, double* results) {
    for (size_t i = 0; i+3 < NR; i+=4) {
        __m256d M_intr[4];
        for (size_t j = 0; j < 4; j++) { M_intr[j] = _mm256_load_pd(Ms + 4*i + 4*j); }

        __m256d v_intr[2];
        for (size_t j = 0; j < 2; j++) { v_intr[j] = _mm256_load_pd(Vs + 2*i + 4*j); }

        __m256d avx_results_intr = mm_vT_M_v_avx2(M_intr[0], M_intr[1],
                                                    M_intr[2], M_intr[3],
                                                    v_intr[0], v_intr[1]);
        _mm256_store_pd(results + i, avx_results_intr);
    }
}
#endif

int main() {

    double* results = static_cast<double*>(aligned_alloc(32, NR * sizeof(double)));
    double* Ms = static_cast<double*>(aligned_alloc(32, NR * 4 * sizeof(double)));
    double* Vs = static_cast<double*>(aligned_alloc(32, NR * 2 * sizeof(double)));

    Benchmark<decltype(&vTMv_base)> bench("v.T @ M T v");

    bench.data_loader = &fill;

    bench.add_function(&vTMv_base, "vTMv base", 9*NR);
#ifdef __AVX2__
    bench.add_function(&vTMv_avx2, "vTMv avx2", 9*NR);
#endif

    // Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(Ms, Vs, results);


    free(results);
    free(Ms);
    free(Vs);
}