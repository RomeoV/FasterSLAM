#include "rdtsc_benchmark.h"
#include <iostream>
#include "add_control_noise.h"

#include "ut.hpp"
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

int main() {
    // Note:
    // add_control_noise and add_control_noise_base are the same at the moment

    // Initialize Input
    double V = 3.0;
    double G = -0.008726646250000001;
    double Q[4] = {0.089999999999999997, 0, 0, 0.0027415567718150069};
    int addnoise = 1;
    double VnGn_base[2];
    double VnGn[2];

    // sanity check
    srand(0);
    add_control_noise_base(V,G,Q,addnoise,VnGn_base);
    srand(0);
    add_control_noise(V,G,Q,addnoise,VnGn);

    expect(that % fabs(VnGn[0] - VnGn_base[0]) < 1.0e-10) << "Speed+Noise";
    expect(that % fabs(VnGn[1] - VnGn_base[1]) < 1.0e-10) << "SteeringAngle+Noise";

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&add_control_noise)> bench("add_control_noise Benchmark");

    double work = 500.0; // best-case in flops

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&add_control_noise_base, "add_control_noise_base", work);
    bench.add_function(&add_control_noise, "add_control_noise", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(V,G,Q,addnoise,VnGn);

    return 0;
}
