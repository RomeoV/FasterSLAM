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


using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`


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

int main() {
    const int N = 1000;
    int* recipient = static_cast<int *>(aligned_alloc(32, N * sizeof(int)));

    srand(0);
    fast_srand(0);
    lehmer64_srand(1);
    wyrng_srand(0);
    wyhash64_srand(0);
    xorshift128plus_srand(1,1);
    pcg32_srand(1,1);
    avx_xorshift128plus_init(1,1);

    auto rand_lambda = random_lambda(&rand);
    auto fastrand_lambda = random_lambda(&fast_rand);
    auto lehmer64_lambda = random_lambda(&lehmer64);
    auto wyrng_lambda = random_lambda(&wyrng);
    auto wyhash64_lambda = random_lambda(&wyhash64);
    auto xorshift128plus_lambda = random_lambda(&xorshift128plus);
    auto pcg32_lambda = random_lambda(&pcg32);
    auto read_rand_lambda = random_lambda(&read_rand);

    //SIMD
    auto avx_xorshift128plus_lambda = random_vec_lambda(&avx_xorshift128plus);
    
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
    bench.add_function(avx_xorshift128plus_lambda, "avx_xorshift128plus", work);
    bench.add_function(read_rand_lambda, "read_rand", work);

    /*
    lehmer64_srand(1);
    fast_srand(0);
    srand(0);
    for ( int i = 0; i<N; i++) {
        std::cout<<rand()<<std::endl;
        std::cout<<fast_rand()<<std::endl;
        std::cout<<lehmer64()<<std::endl<<std::endl;
    }
    */
    

    bench.run_benchmark(recipient, N);

}