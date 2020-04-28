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
    "functional equality"_test = [&] {
        auto is_close = [](auto lhs, auto rhs) -> bool {return lhs - rhs < 1e-14 or rhs - lhs < 1e-14;};
        for (size_t i = 0; i < N; i++) {
            double method1 = pi_to_pi(angles[i]);
            double method2 = pi_to_pi_fmod(angles[i]);
            expect(is_close(method1, method2));
        }
    };

    double(*lambda_pi_to_pi)(double* angles) = [](double* angles){
        double sum=0;
        for (int i = 0; i<N;i++) {
            sum+=pi_to_pi(angles[i]);
        }
        return sum;
    };

    double(*lambda_pi_to_pi_base)(double* angles) = [](double* angles){
        double sum=0;
        for (int i = 0; i<N;i++) {
            sum+=pi_to_pi_base(angles[i]);
        }
        return sum;
    };

    double(*lambda_pi_to_pi_fmod)(double* angles) = [](double* angles){
        double sum = 0;
        for (int i = 0; i<N;i++) {
            sum+=pi_to_pi_fmod(angles[i]);
        }
        return sum;
    };

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&pi_to_pi)> bench("pi_to_pi Benchmark");

    double work = 5.0; // best-case

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&pi_to_pi_base, "pi_to_pi_base", work);
    bench.add_function(&pi_to_pi, "pi_to_pi", work);
    bench.add_function(&pi_to_pi_fmod, "pi_to_pi_fmod", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(angles[0]);

    
    Benchmark<decltype(lambda_pi_to_pi)> bench_lambda("pi_to_pi on array Benchmark");

    // Add lambda functions to aggregate over range of inputs.
    bench_lambda.add_function(lambda_pi_to_pi_base, "pi_to_pi_base", work*N);
    bench_lambda.add_function(lambda_pi_to_pi, "pi_to_pi", work);
    bench_lambda.add_function(lambda_pi_to_pi_fmod, "pi_to_pi_fmod", work*N);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench_lambda.run_benchmark(angles);

    //bench_lambda.summary_long();

    
    /*
    // Alternative (much slower here, but nicer to look at. Generally useful if you want to average over a few inputs). Yields averages over all runs.
    
    Benchmark<decltype(&pi_to_pi)> multi_bench("pi_to_pi Benchmark");

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    multi_bench.add_function(&pi_to_pi, "pi_to_pi", 6);
    multi_bench.add_function(&pi_to_pi_fmod, "pi_to_pi_fmod", 6);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    for (int i = 0; i<N; i++) {
        // You could set the data_loader function here to generate new input. multi_bench.data_loader =&my_load_func_i...
        multi_bench.run_benchmark(angles[i]);
    }
    */

    return 0;
}