#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "get_observations.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

auto data_loader(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, int *idf, size_t *nidf, Vector2d z[]) {
    for (int i = 0; i < *nidf; i++) {
        idf[i] = i;
    }
}

int main() {
    cVector3d x = {0.674090417131751, -0.030904471130924017, -0.0073589032333721818};
    const double rmax = 30; 
    const size_t lm_rows = 35;
    size_t exact_nidf = 35, nidf = 35;
    double lm[70];
    int exact_idf[35], idf[35];
    Vector2d exact_z[2], z[2];

    // Fill lm from file
    FILE* fp = fopen("inputfiles_test/lm.txt", "r");  expect(fp != 0);
    for (size_t i = 0; i < lm_rows; i++) {
        fscanf(fp, "%lf\t%lf\n", &lm[i*2+0], &lm[i*2+1]);
    }
    fclose(fp);

    // Test:  
    data_loader(x, rmax, lm, lm_rows, exact_idf, &exact_nidf, exact_z);
    get_observations_base(x, rmax, lm, lm_rows, exact_idf, &exact_nidf, exact_z);
    
    data_loader(x, rmax, lm, lm_rows, idf, &nidf, z);
    get_observations(x, rmax, lm, lm_rows, idf, &nidf, z);
    
    // Check x
    expect( exact_nidf == nidf );
    for (int i = 0; i < 2; i++) {
        double error0 = fabs( z[i][0] - exact_z[i][0] );
        double error1 = fabs( z[i][1] - exact_z[i][1] );
        expect(that % error0 < 1e-12) << i;
        expect(that % error1 < 1e-12) << i;
    }

    Benchmark<decltype(&get_observations)> bench("get_observations benchmark");
    double work = 15*lm_rows; // TODO Count work //best-case
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&get_observations_base, "base", work);
    bench.add_function(&get_observations, "active", work);

    bench.run_benchmark(x, rmax, lm, lm_rows, idf, &nidf, z);

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

