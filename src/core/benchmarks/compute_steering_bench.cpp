#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "compute_steering.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`


auto data_loader(cVector3d x, double* wp, const size_t N_wp, const double minD, 
                      const double rateG, const double maxG, const double dt, 
                      int* iwp, double* G) {
    *iwp=1;
    *G = 0.0;
}
int main() {
    double x[3] = {0,0,0};
    const int N_wp = 3;
    double wp[2*N_wp] = {0,0,1,1,2,2};
    double minD = 12.0;
    double rateG = 0.1;
    double maxG = 1.5;
    double dt = 1.0;
    int iwp = 1;
    double G = 0.0;

    int iwp_exact = 1;
    double G_exact = 0.0;

    data_loader(x, wp, N_wp, minD, rateG, maxG, dt, &iwp_exact, &G_exact);
    compute_steering_base(x, wp, N_wp, minD, rateG, maxG, dt, &iwp_exact, &G_exact);
    data_loader(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);
    compute_steering(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);

    expect(that % fabs(G - G_exact) < 1.0e-10) << "G";
    expect(that % (iwp - iwp_exact) == 0 ) << "iwp";

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
//    Benchmark<decltype(&compute_steering_base)> bench("compute_steering Benchmark");
//
//    data_loader(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);
//    bench.data_loader = data_loader;
//    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
//    // First function should always be the base case you want to benchmark against!
//    bench.add_function(&compute_steering_base, "base", 0.0);
//    bench.funcFlops[0] = compute_steering_base_flops(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);
//    bench.funcBytes[0] = 8*compute_steering_base_memory(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);
//    bench.add_function(&compute_steering, "active", 0.0);
//    bench.funcFlops[1] = compute_steering_active_flops(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);
//    bench.funcBytes[1] = 8*compute_steering_active_memory(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);
//
//    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
//    bench.run_benchmark(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);

    return 0;
}
