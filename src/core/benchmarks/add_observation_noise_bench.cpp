#include "rdtsc_benchmark.h"

/* Jonathan's version
#include <iostream>
#include "add_observation_noise.h"

int main() {
    // Initialize Input
    Vector2d z[2] = {{25.77444969441645, -1.4733774573801159}, {25.276107769232738, 0.13836675004531551}};
    const int zlen = 2;
    double R[4] = {0.010000000000000002, 0, 0, 0.00030461741909055634};
    const int addnoise = 1;

    const size_t N = 1000;

    double(*lambda)(Vector2d* z, int zlen, double *R, int addnoise) = [](Vector2d* z, int zlen, double *R, int addnoise){
        double sum=0;
        for (int i = 0; i<N; i++) {
            add_observation_noise(z, zlen, R, addnoise);
        }
        return sum;
    };

    double(*lambda_base)(Vector2d* z, int zlen, double *R, int addnoise) = [](Vector2d* z, int zlen, double *R, int addnoise){
        double sum=0;
        for (int i = 0; i<N;i++) {
            add_observation_noise(z, zlen, R, addnoise);
        }
        return sum;
    };

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&add_observation_noise)> bench("add_observation_noise Benchmark");

    double work = 50000.0; // best-case in flops

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&add_observation_noise, "add_observation_noise", work);
    bench.add_function(&add_observation_noise_base, "add_observation_noise_base", work);
    //bench.add_function(&add_observation_noise_fmod, "add_observation_noise_fmod", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(z, zlen, R, addnoise);

    
    Benchmark<decltype(lambda)> bench_lambda("add_observation_noise on array Benchmark");

    // Add lambda functions to aggregate over range of inputs.
    bench_lambda.add_function(lambda, "add_observation_noise", work*N);
    bench_lambda.add_function(lambda_base, "aadd_observation_noise_base", work*N);
    //bench_lambda.add_function(lambda_fmod, "add_observation_noise_fmod", work*N);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench_lambda.run_benchmark(z, zlen, R, addnoise);


    // Free memory
    // destroy(A);
    // destroy(x);
    // destroy(y);

    return 0;
}
*/

#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "add_observation_noise.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

auto data_loader(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise) {
    Vector2d z_init[2] = {{25.77444969441645, -1.4733774573801159}, {25.276107769232738, 0.13836675004531551}};
    //std::copy(z_init[0], z_init[0]+2, z[0]);
    //std::copy(z_init[1], z_init[1]+2, z[1]);
    std::copy(&z_init[0][0], &z_init[0][0]+4, &z[0][0]);
    // reset seed
    srand(1994);
}

int main() {

    // Test:
    Vector2d exact_z[2], z[2];
    const int zlen = 2;
    cMatrix2d R = {0.010000000000000002, 0, 0, 0.00030461741909055634};
    const int addnoise = 1;

    data_loader(exact_z, zlen, R, addnoise);
    add_observation_noise_base(exact_z, zlen, R, addnoise);
    
    data_loader(z, zlen, R, addnoise);
    add_observation_noise(z, zlen, R, addnoise);
    
    // Check x
    for (int i = 0; i < zlen; i++) {
        double error0 = fabs( z[i][0] - exact_z[i][0] );
        double error1 = fabs( z[i][1] - exact_z[i][1] );
        expect(that % error0 < 1e-12) << i;
        expect(that % error1 < 1e-12) << i;
    }

    Benchmark<decltype(&add_observation_noise)> bench("add_observation_noise benchmark");
    double work = 6*zlen +2; // TODO Count work
    bench.data_loader = data_loader; // To guarantee same inputs
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&add_observation_noise_base, "base", work);
    bench.add_function(&add_observation_noise, "active", work);

    bench.run_benchmark(z, zlen, R, addnoise);

    return 0;
}