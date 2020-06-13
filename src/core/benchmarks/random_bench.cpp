#include "rdtsc_benchmark.h"
#include "ut.hpp"
#include <immintrin.h>
#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>
#include <random>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "fastrand.h"
#include "linalg.h"


using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

/*

O3:
 |  i | Name                           |      Work [flops] |        Time [cyc] | Perf. [flops/cyc] |       Speedup [-] | 
 |  0 | random_base                    |              1000 |        27306.0762 |            0.0366 |            1.0000 | 
 |  1 | fastrand                       |              1000 |         7033.3294 |            0.1422 |            3.8824 | 
 |  2 | lehmer64                       |              1000 |         7069.1736 |            0.1415 |            3.8627 | 
 |  3 | wyrng                          |              1000 |         6185.3129 |            0.1617 |            4.4147 | 
 |  4 | wyhash64                       |              1000 |         7079.0271 |            0.1413 |            3.8573 | 
 |  5 | pcg32                          |              1000 |         8458.2665 |            0.1182 |            3.2283 | 
 |  6 | xorshift128plus                |              1000 |         7226.9642 |            0.1384 |            3.7784 | 
 |  7 | avx_xorshift128plus            |              1000 |         2826.4394 |            0.3538 |            9.6609 | 
 |  8 | xorshf96                       |              1000 |         6794.7543 |            0.1472 |            4.0187 | 
 |  9 | read_rand                      |              1000 |         8410.5777 |            0.1189 |            3.2466 | 
 
O0:
fastrand Benchmark (1 runs, avg):
 |  i | Name                           |      Work [flops] |        Time [cyc] | Perf. [flops/cyc] |       Speedup [-] | 
 |  0 | random_base                    |              1000 |        27911.7286 |            0.0358 |            1.0000 | 
 |  1 | fastrand                       |              1000 |         9373.5196 |            0.1067 |            2.9777 | 
 |  2 | lehmer64                       |              1000 |        10981.6625 |            0.0911 |            2.5417 | 
 |  3 | wyrng                          |              1000 |        22765.0730 |            0.0439 |            1.2261 | 
 |  4 | wyhash64                       |              1000 |        22940.7960 |            0.0436 |            1.2167 | 
 |  5 | pcg32                          |              1000 |        16455.2145 |            0.0608 |            1.6962 | 
 |  6 | xorshift128plus                |              1000 |        16270.4372 |            0.0615 |            1.7155 | 
 |  7 | avx_xorshift128plus            |              1000 |        13071.6596 |            0.0765 |            2.1353 | 
 |  8 | xorshf96                       |              1000 |        19125.6207 |            0.0523 |            1.4594 | 
 |  9 | read_rand                      |              1000 |         9839.0255 |            0.1016 |            2.8368 | 
 
 
 */

template<typename func>
std::function<void (int *, int)> random_lambda(func rand_func) {
    return [=] (int* recipient, int N)
    {
        for (int i = 0; i<N; i++) {
            recipient[i] = rand_func();
            
        }
    };
}



template<typename func>
std::function<void (int *, int)> random_vec_lambda(func rand_func) {
    auto mask = _mm256_set1_epi32 (1);
    return [=] (int* recipient, int N)
    {
        for (int i = 0; i<N; i+=4) {
            _mm256_maskstore_epi32(recipient+i,mask, rand_func());
        }
    };
}


template<typename func>
std::function<void (int *, int)> random_vec32_lambda(func rand_func) {
    auto mask = _mm256_set1_epi32 (1);
    return [=] (int* recipient, int N)
    {
        for (int i = 0; i<N; i+=8) {
            _mm256_maskstore_epi32(recipient+i,mask, rand_func());
        }
    };
}


int main() {
#ifdef __AVX2__
    const int N = 10000;
    int* recipient = static_cast<int *>(aligned_alloc(32, N * sizeof(int)));

    srand(0);
    fast_srand(0);
    lehmer64_srand(1);
    wyrng_srand(0);
    wyhash64_srand(0);
    xorshift128plus_srand(1,1);
    pcg32_srand(1,1);
    avx_xorshift128plus_init(1,1);

    avx_fast_srand(0,10,100,1000);

    uint64_t init_state[8] __attribute__((aligned(32))) = {1,1,1,1,1,1,1,1};
    uint64_t init_seq[8] __attribute__((aligned(32))) = {1,2,3,4,5,6,7,8};


    avx2_pcg32_srand(init_state, init_seq);


    auto rand_lambda = random_lambda(&rand);
    auto fastrand_lambda = random_lambda(&fast_rand);
    auto lehmer64_lambda = random_lambda(&lehmer64);
    auto wyrng_lambda = random_lambda(&wyrng);
    auto wyhash64_lambda = random_lambda(&wyhash64);
    auto xorshift128plus_lambda = random_lambda(&xorshift128plus);
    auto pcg32_lambda = random_lambda(&pcg32);
    auto xorshf96_lambda = random_lambda(&xorshf96);
    auto read_rand_lambda = random_lambda(&read_rand);

    //SIMD
    auto avx_xorshift128plus_lambda = random_vec_lambda(&avx_xorshift128plus);
    auto avx_fast_rand_lambda = random_vec_lambda(&avx_fast_rand);
    auto avx2_pcg32_lambda = random_vec32_lambda(&avx2_pcg32);

    Benchmark<decltype(fastrand_lambda)> bench("fastrand Benchmark");

    double work = N; // best-case in flops

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(rand_lambda, "random_base", work);
    bench.add_function(fastrand_lambda, "fastrand", work);
    bench.add_function(lehmer64_lambda, "lehmer64", work);
    bench.add_function(wyrng_lambda, "wyrng", work);
    bench.add_function(wyhash64_lambda, "wyhash64", work);
    bench.add_function(pcg32_lambda, "pcg32", work);
    bench.add_function(xorshift128plus_lambda, "xorshift128plus", work);
    
    bench.add_function(xorshf96_lambda, "xorshf96", work);
    bench.add_function(read_rand_lambda, "read_rand", work);

    bench.add_function(avx_xorshift128plus_lambda, "avx_xorshift128plus", work);
    bench.add_function(avx_fast_rand_lambda, "avx_fast_rand", work);
    bench.add_function(avx2_pcg32_lambda, "avx2_pcg32", work);
    bench.run_benchmark(recipient, N);

    free(recipient);

#else
#warning "Disabled trigonometry_bench because AVX2 is not supported!"
#endif
    return 0;
}
