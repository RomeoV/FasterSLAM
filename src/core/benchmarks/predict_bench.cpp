#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>
#include "particle.h"
#include "configfile.h"

#include "predict.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

void data_loader(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt) {
    particle->xv[0] = 1;
    particle->xv[1] = 1;
    particle->xv[2] = M_PI/2;
}

int main() {
    SWITCH_PREDICT_NOISE = 0;
    double xv_initial[3] = {0.0,0.0,0.0};
    double Q[4] = {0.01,0,0,0.03};

    double V = 1.0;
    double G = 2.0;
    double dt = 0.1;
    double WB = 0.1;
    
    Particle p;
    initParticle(&p, 1, xv_initial);
    data_loader(&p, V, G, Q, WB, dt);

    Particle p_exact;
    initParticle(&p_exact, 1, xv_initial);
    data_loader(&p_exact, V, G, Q, WB, dt);

    predict(&p, V, G, Q, WB, dt),
    predict_base(&p_exact, V, G, Q, WB, dt);
    for (int i = 0; i < 3; i++) {
        expect(that % fabs(p.xv[i]-p_exact.xv[i])<= 1.0e-10) << i;
    }

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&predict_base)> bench("predict Benchmark");    
    
    data_loader(&p, V, G, Q, WB, dt); // might be enough to call it above
    bench.data_loader = data_loader;
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&predict_base, "base", 0.0);
    bench.funcFlops[0] = predict_base_flops(&p, V, G, Q, WB, dt);
    bench.funcBytes[0] = predict_base_memory(&p, V, G, Q, WB, dt);

    // predict simply calls predict_base at the moment
    bench.add_function(&predict, "active", 0.0);
    bench.funcFlops[1] = predict_base_flops(&p, V, G, Q, WB, dt);
    bench.funcBytes[1] = predict_base_memory(&p, V, G, Q, WB, dt);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(&p, V, G, Q, WB, dt);

    return 0;
}