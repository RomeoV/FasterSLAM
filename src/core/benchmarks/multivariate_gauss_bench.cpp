#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "multivariate_gauss.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

auto data_loader(cVector2d x, cMatrix2d P, Vector2d result) {
    // empty
}

int main() {

    // Test: 
    cVector2d x = {3.0, -0.008726646250000001};
    cMatrix2d P = {0.089999999999999997, 0, 0, 0.0027415567718150069};
    Vector2d exact_result, result;
    
    srand(1994);
    data_loader(x, P, exact_result);
    multivariate_gauss(x, P, exact_result);
    
    srand(1994);
    data_loader(x, P, result);
    multivariate_gauss(x, P, result);
    
    // Check x
    for (int i = 0; i < 2; i++) {
        double error = fabs( result[i] - exact_result[i] );
        expect(that % error < 1e-12) << i;
    }

    Benchmark<decltype(&multivariate_gauss)> bench("multivariate_gauss benchmark");
    double work = 17;
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&multivariate_gauss_base, "base", work);
    bench.add_function(&multivariate_gauss, "active", work);

    bench.run_benchmark(x, P, result);

    /*
    // Alternative (much slower here, but nicer to look at. Generally useful if you want to average over a few inputs). Yields averages over all runs.
    
    Benchmark<decltype(&pi_to_pi)> bench("pi_to_pi Benchmark");

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&pi_to_pi, "pi_to_pi", 6);
    bench.add_function(&pi_to_pi_fmod, "pi_to_pi_fmod", 6);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    for (int i = 0; i<N; i++) {
        // You could set the data_loader function here to generate new input. bench.data_loader =&my_load_func_i...
        bench.run_benchmark(angles[i]);
    }
    */

    return 0;
}

