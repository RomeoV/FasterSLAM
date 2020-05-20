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


int main() {

    // Test: 
    cVector2d x = {3.0, -0.008726646250000001};
    cMatrix2d P = {0.089999999999999997, 0, 0, 0.0027415567718150069};
    Vector2d exact_result, result;
    
    srand(1994);
    multivariate_gauss_base(x, P, exact_result);
    
    srand(1994);
    multivariate_gauss(x, P, result);
    
    // Check x
    for (int i = 0; i < 2; i++) {
        double error = fabs( result[i] - exact_result[i] );
        expect(that % error < 1e-12) << i;
    }

    Benchmark<decltype(&multivariate_gauss)> bench("multivariate_gauss benchmark");
    double work = 17;
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&multivariate_gauss_base, "base", work);
    //bench.add_function(&multivariate_gauss_fast_rand, "fast pseudo RNG", work);
    bench.add_function(&multivariate_gauss, "active", work);

    bench.run_benchmark(x, P, result);

    return 0;
}

