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


using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`


void fill_angles(double* angles, int N, double lower_bound, double upper_bound) {
    for (int i = 0; i<N; i++) {
        angles[i] = lower_bound + (upper_bound -lower_bound) / (N-1) * i;
    }
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

int main() {
    init_sin();
    init_sin2();
    const int N = 10000;
    double angles[N];
    double vec_angles[N];
    double lower_bound = 0.0;
    double upper_bound = M_PI;

    fill_angles(angles, N, lower_bound, upper_bound);
    fill_angles(vec_angles, N, 0, 0); // For some reason this is required

    // Compare two methods
    "functional equality"_test = [&] {
        auto is_close = [&](auto lhs, auto rhs) -> bool {return lhs - rhs < 1e-8 or rhs - lhs < 1e-8;};
        
        for (int i = 0; i<N; i+=4) {
            auto vec = _mm256_load_pd(angles+i);
            _mm256_store_pd(vec_angles+i, read_sin_vec(vec));
        }

        
        for (int i = 0; i < N; i++) {
            double method1 = std::sin(angles[i]);
            double method2 = read_sin(angles[i]);
            double method3 = tscheb_sine(angles[i]);
            expect(is_close(method1, method2));
            expect(is_close(method1, method3));
            expect(is_close(method1, vec_angles[i]));
        }
    };
    
    auto std_sin_lambda = trig_lambda(&sine);
    auto read_sin_lambda = trig_lambda(&read_sin);
    auto tscheb_sin_lambda = trig_lambda(&tscheb_sine);

    auto read_sin_vec_lambda = trig_vec_lambda(&read_sin_vec);
    auto read_sin2_vec_lambda = trig_vec_lambda(&read_sin2_vec);


    Benchmark<decltype(std_sin_lambda)> bench("Trigonometry Benchmark");

    // Add lambda functions to aggregate over range of inputs.
    bench.add_function(std_sin_lambda, "base", N);
    bench.add_function(read_sin_lambda, "read_sin", N);
    bench.add_function(tscheb_sin_lambda, "tscheb_sine",N);

    bench.add_function(read_sin_vec_lambda, "read_sin_vec", N);
    bench.add_function(read_sin2_vec_lambda, "read_sin2_vec", N);

    auto loader = data_loader_lambda(fill_angles, lower_bound, upper_bound);

    bench.data_loader = loader;

    bench.run_benchmark(angles, N);





}
