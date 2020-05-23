#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "trigonometry.h"
#include "linalg.h"
#include <immintrin.h>
#include "tscheb_sine.h"
#include "typedefs.h"


using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`



void fill_angles(double* angles, int N, double lower_bound, double upper_bound) {
    // for (int i = 0; i<N; i++) {
    //     angles[i] = lower_bound + (upper_bound -lower_bound) / (N-1) * i;
    // }
    srand(0);
    fill_rand(angles, N, lower_bound, upper_bound);
}

template<typename func>
std::function<void (double*, int)> data_loader_lambda(func trig_func, double lower_bound, double upper_bound) {
    return [=] (double* angles, int N)
    {
        fill_angles(angles, N, lower_bound, upper_bound);
    };
}



template<typename func>
std::function<void (double*, int)> trig_lambda(func trig_func) {
    return [=] (double* angles, int N)
    {
        for (int i = 0; i<N; i++) {
            angles[i] = trig_func(angles[i]);
        }
    };
}

template<typename func>
std::function<void (double*, int)> trig_vec_lambda(func trig_func) {
    return [=] (double* angles, int N)
    {
        for (int i = 0; i<N; i+=4) {
            auto vec = _mm256_load_pd(angles+i);
            _mm256_store_pd(angles+i, trig_func(vec));
        }
    };
}



inline double sine(double angle) {
    return std::sin(angle);
}

inline double tscheb_sine(double angle) {
    return tscheb_dsine(angle, true);
}

inline double tscheb_sine_nn(double angle) {
    return tscheb_dsine(angle, false);
}
int main() {
#ifdef __AVX2__
    
    const int N = 100000;
    double angles[N]  __attribute__((aligned(32)));
    double vec_angles[N]  __attribute__((aligned(32)));
    double lower_bound = -3* M_PI;
    double upper_bound = 3* M_PI;

    init_sin();
    init_sin2();

    fill_angles(angles, N, lower_bound, upper_bound);
    fill_angles(vec_angles, N, 0, 0); // For some reason this is required

    // Compare two methods
    "functional equality"_test = [&] {
        auto is_close = [&](auto lhs, auto rhs) -> bool {return std::abs(lhs - rhs) < 1e-6;};
        
        for (int i = 0; i<N; i+=4) {
            auto vec = _mm256_load_pd(angles+i);
            _mm256_store_pd(vec_angles+i, tscheb_sin_avx(vec));
        }

        
        for (int i = 0; i < N; i++) {
            double method1 = std::sin(angles[i]);
            //double method2 = read_sin(angles[i]);
            double method3 = tscheb_sin(angles[i]);
            //expect(is_close(method1, method2))<< "read_sin" << method1 << method2;
            expect(is_close(method1, method3))<< "tscheb";
            expect(is_close(method1, vec_angles[i]))<< "read_sin_vec"<< method1 << vec_angles[i]<< i << angles[i];
        }
    };
    
    auto std_sin_lambda = trig_lambda(&sine);
    auto read_sin_lambda = trig_lambda(&read_sin);
    auto tscheb_sin_lambda = trig_lambda(&tscheb_sine);
    auto tscheb_sin_simple_lambda = trig_lambda(&tscheb_sin);
    auto tscheb_sin_not_normalized_lambda = trig_lambda(&tscheb_sine_nn);

    auto read_sin_vec_lambda = trig_vec_lambda(&read_sin_vec);
    auto read_sin2_vec_lambda = trig_vec_lambda(&read_sin2_vec);
    auto tscheb_sin_vec_lambda = trig_vec_lambda(&tscheb_sin_avx);


    Benchmark<decltype(std_sin_lambda)> bench("Trigonometry Benchmark");

    // Add lambda functions to aggregate over range of inputs. 
    // there is no instrumentation available
    bench.add_function(std_sin_lambda, "base", N);
    bench.funcFlops[0] = N * tp.sin;
    bench.funcBytes[0] = N;
    bench.add_function(read_sin_lambda, "read_sin", N);
    bench.funcFlops[1] = 1;
    bench.funcBytes[1] = 1;

    bench.add_function(tscheb_sin_lambda, "tscheb_sine_normalized",N);
    bench.funcFlops[2] = N * tscheb_dsine_flops(3, true);
    bench.funcBytes[2] = N;

     bench.add_function(tscheb_sin_simple_lambda , "tscheb_sine_partially_normalized",N);
    bench.funcFlops[3] = N * tscheb_dsine_flops(3, true)+ 2* N * tp.add;
    bench.funcBytes[3] = N;

    bench.add_function(tscheb_sin_not_normalized_lambda, "tscheb_sine_not_normalized", N);
    bench.funcFlops[4] = N * tscheb_dsine_flops(4, false);
    bench.funcBytes[4] = N;

    bench.add_function(read_sin_vec_lambda, "read_sin_avx", N);
    bench.funcFlops[5] = 1;
    bench.funcBytes[5] = N;
    bench.add_function(read_sin2_vec_lambda, "read_sin_symmetry_avx", N);
    bench.funcFlops[6] = 1;
    bench.funcBytes[6] = N;
    bench.add_function(tscheb_sin_vec_lambda, "tscheb_sin_avx_partially_normalized",N);
    bench.funcFlops[7] = N/4 * tscheb_sin_avx_flops(_mm256_load_pd(angles));
    bench.funcBytes[7] = N;

    auto loader = data_loader_lambda(fill_angles, lower_bound, upper_bound);

    bench.data_loader = loader;

    bench.run_benchmark(angles, N);
#else
#warning "Disabled trigonometry_bench because AVX2 is not supported!"
#endif
    return 0;
}

