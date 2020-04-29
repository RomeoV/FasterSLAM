#include "rdtsc_benchmark.h"
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

