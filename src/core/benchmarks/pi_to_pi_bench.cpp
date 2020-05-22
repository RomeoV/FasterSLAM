#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "pi_to_pi.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

template <typename F>
struct NamedFunction {
    std::string name;
    F function;
};

int main() {
    // Initialize Input
    const size_t N = 1000;
    double angles[N];
    std::generate(angles, angles+N, []{return 2*M_PI*((rand()%1000)/100. - 5);});

    // Compare two methods
    "functional equality base-fmod"_test = [&] {
        auto is_close = [](auto lhs, auto rhs) -> bool {return lhs - rhs < 1e-14 or rhs - lhs < 1e-14;};
        for (size_t i = 0; i < N; i++) {
            double method1 = pi_to_pi_base(angles[i]);
            double method2 = pi_to_pi_fmod(angles[i]);
            expect(is_close(method1, method2));
        }
    };

    "functional equality base-opt"_test = [&] {
        auto is_close = [](auto lhs, auto rhs) -> bool {return lhs - rhs < 1e-14 or rhs - lhs < 1e-14;};
        for (size_t i = 0; i < N; i++) {
            double method1 = pi_to_pi_base(angles[i]);
            double method2 = pi_to_pi(angles[i]);
            expect(is_close(method1, method2));
        }
    };

    double(*lambda_pi_to_pi)(double* angles) = [](double* angles){
        double sum=0;
        for (int i = 0; i<N; i++) {
            sum+=pi_to_pi(angles[i]);
        }
        return sum;
    };

    double(*lambda_pi_to_pi_base)(double* angles) = [](double* angles){
        double sum=0;
        for (int i = 0; i<N; i++) {
            sum+=pi_to_pi_base(angles[i]);
        }
        return sum;
    };

    double(*lambda_pi_to_pi_fmod)(double* angles) = [](double* angles){
        double sum = 0;
        for (int i = 0; i<N; i++) {
            sum+=pi_to_pi_fmod(angles[i]);
        }
        return sum;
    };

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&pi_to_pi)> bench("pi_to_pi Benchmark");

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&pi_to_pi_base, "pi_to_pi_base", 0.0);
    bench.funcFlops[0] = pi_to_pi_base_flops(angles[0]);
    bench.funcBytes[0] = pi_to_pi_base_memory(angles[0]);
    bench.add_function(&pi_to_pi, "pi_to_pi", 0.0);
    bench.funcFlops[1] = pi_to_pi_active_flops(angles[0]);
    bench.funcBytes[1] = pi_to_pi_active_memory(angles[0]);
    bench.add_function(&pi_to_pi_fmod, "pi_to_pi_fmod", 0.0);
    bench.funcFlops[2] = pi_to_pi_active_flops(angles[0]);
    bench.funcBytes[2] = pi_to_pi_active_memory(angles[0]);
    bench.add_function(&pi_to_pi_nongeneral, "pi_to_pi_nongeneral", 0.0); // 2
    bench.funcFlops[3] = pi_to_pi_active_flops(angles[0]);
    bench.funcBytes[3] = pi_to_pi_active_memory(angles[0]);
    bench.add_function(&pi_to_pi_while, "pi_to_pi_while", 0.0); // 4 // on average
    bench.funcFlops[4] = pi_to_pi_active_flops(angles[0]);
    bench.funcBytes[4] = pi_to_pi_active_memory(angles[0]);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(angles[0]);

    
    Benchmark<decltype(lambda_pi_to_pi)> bench_lambda("pi_to_pi on array Benchmark");

    // Add lambda functions to aggregate over range of inputs.
    double work = 0.0;
    double memory = 0.0;
    bench_lambda.add_function(lambda_pi_to_pi_base, "pi_to_pi_base", 0.0);
    for(int i=0; i<N; i++){
        work += pi_to_pi_base_flops(angles[i]);
        memory += pi_to_pi_base_memory(angles[i]);
    }
    bench.funcFlops[0] = work;
    bench.funcBytes[0] = memory;
    bench_lambda.add_function(lambda_pi_to_pi, "pi_to_pi", 0.0);
    for(int i=0; i<N; i++){
        work += pi_to_pi_active_flops(angles[i]);
        memory += pi_to_pi_active_memory(angles[i]);
    }
    bench.funcFlops[1] = work;
    bench.funcBytes[1] = memory;
    bench_lambda.add_function(lambda_pi_to_pi_fmod, "pi_to_pi_fmod", 0.0);
     for(int i=0; i<N; i++){
        work += pi_to_pi_active_flops(angles[i]);
        memory += pi_to_pi_active_memory(angles[i]);
    }
    bench.funcFlops[2] = work;
    bench.funcBytes[2] = memory;

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench_lambda.run_benchmark(angles);

    return 0;
}