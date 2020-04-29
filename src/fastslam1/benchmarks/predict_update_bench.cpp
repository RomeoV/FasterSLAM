#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "particle.h"
#include "fastslam1_utils.h"
#include "configfile.h"

#include "predict_update.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`


void data_loader(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    srand(0);
    xtrue[0]= 0.0;
    xtrue[1]= 0.0;
    xtrue[2]= 0.0;
    *iwp=1;
    *G = 0.0;

    for (int i = 0; i< N; i++) {
        particles[i].xv[0] = xtrue[0];
        particles[i].xv[1] = xtrue[1];
        particles[i].xv[2] = xtrue[2];
    }
}
int main() {
    const int N = 100;
    double minD = 0.1;
    double rateG = 1.0;
    double maxG = 1.5;
    double dt = 1.0;
    double V = 3.0;
    setup_initial_Q_R();  // Remove this soon!!!

    int N_waypoints = 3;
    double wp[6] = {0,0,1,1,2,2};

    Particle particles[N];
    int iwp = 1;
    double G = 0.0;
    Vector3d xtrue;
    data_loader(wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles);
    predict_update(wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G, particles);

    Particle particles_exact[N];
    int iwp_exact = 1;
    double G_exact = 0.0;
    Vector3d xtrue_exact;
    data_loader(wp, N_waypoints, V, *Q, dt, N, xtrue_exact, &iwp_exact, &G_exact, particles_exact);
    predict_update_base(wp, N_waypoints, V, *Q, dt, N, xtrue_exact, &iwp_exact, &G_exact, particles_exact);
    

    expect(that % fabs(G - G_exact) < 1.0e-10) << "G";
    expect(that % (iwp - iwp_exact) == 0 ) << "iwp";

    for (int i =0; i<N; i++) {
        expect(that % fabs(particles[i].xv[0] - particles_exact[i].xv[0]) < 1.0e-10);
        expect(that % fabs(particles[i].xv[1] - particles_exact[i].xv[1]) < 1.0e-10);
        expect(that % fabs(particles[i].xv[2] - particles_exact[i].xv[2]) < 1.0e-10);
    }

    expect(that % fabs(xtrue[0] - xtrue_exact[0]) < 1.0e-10);
    expect(that % fabs(xtrue[1] - xtrue_exact[1]) < 1.0e-10);
    expect(that % fabs(xtrue[2] - xtrue_exact[2]) < 1.0e-10);

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&predict_update_base)> bench("predict_update Benchmark");

    double work = 41 + N*19; // best

    bench.data_loader = data_loader;
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&predict_update_base, "base", work);
    bench.add_function(&predict_update, "active", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles);

    G= M_PI;
    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles);

    minD = 5.0;
    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles);

    bench.destructor_output = false;
    bench.summary_long();
    return 0;
}