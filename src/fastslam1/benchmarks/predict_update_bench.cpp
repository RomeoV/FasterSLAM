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
#include "linalg.h"
#include "trigonometry.h"
#include "fastrand.h"

#include "predict_update.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`


void data_loader(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    srand(0);
    
    uint64_t init_state[8] = {1,1,1,1,1,1,1,1};
    uint64_t init_seq[8] = {1,2,3,4,5,6,7,8};
    pcg32_srand(1,1);

    avx2_pcg32_srand(init_state, init_seq);
    avx_xorshift128plus_init(1,1);
    xtrue[0]= 1.0;
    xtrue[1]= 2.0;
    xtrue[2]= 0.0;
    *iwp=1;
    *G = 0.0;

    for (int i = 0; i< N; i++) {
        particles[i].xv[0] = xtrue[0];
        particles[i].xv[1] = xtrue[1];
        particles[i].xv[2] = xtrue[2];
    }
}

void init_particles(Particle* particle, double* xv, double* Pv, const size_t N) {
    for (int i = 0; i<N; i++) {
        initParticle(particle+i, 0, xv+3*i);
    } 
}

void init_particles_contigous(Particle* particle, double* xv, double* Pv, const size_t N) {
    for (int i = 0; i<N; i++) {
        particle[i].xv = xv+3*i;
        particle[i].Pv = Pv+9*i;
        initParticle_prealloc(particle+i, 0, xv+3*i);
    } 
}

// I will try to add this as smooth as possible to the benchmark, but for now do this
void set_work(Benchmark<decltype(&predict_update)>& bench, 
                double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    data_loader(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);
    bench.funcFlops[0] = predict_update_base_flops(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);
    data_loader(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);
    bench.funcBytes[0] = 8 * predict_update_base_memory(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);
    data_loader(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);
    for (int i = 1; i < bench.numFuncs; i++) {
        bench.funcFlops[i] = predict_update_active_flops(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);
        data_loader(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);
        bench.funcBytes[i] = 8* predict_update_active_memory(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);
        data_loader(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);
    }
}


void cleanup_members(Particle* particles, int N) {
    for(int i = 0; i<N; i++) {
        delParticleMembers_prealloc(particles+i);
    }
}

int main() {
    //init_sin();
    //init_sin2();
    SWITCH_PREDICT_NOISE=1;
    SWITCH_CONTROL_NOISE=1;
    const int N = 100;
    double minD = 0.1;
    double rateG = 1.0;
    double maxG = 1.5;
    double dt = 1.0;
    double V = 3.0;
    setup_initial_Q_R();  // Remove this soon!!!

    int N_waypoints = 3;
    double wp[6] = {0,0,1,1,2,2};
    double xv[3*N] __attribute__ ((aligned(32)));
    double Pv[9*N] __attribute__ ((aligned(32)));
    fill(xv, 2*N, 0.0);

    Particle particles[N];
    init_particles_contigous(particles, xv, Pv, N);
    int iwp = 1;
    double G = 0.0;
    Vector3d xtrue;
    data_loader(wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles);
    predict_update_fast_normal_rand(wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G, particles);

    
    double xv_exact[3*N] __attribute__ ((aligned(32)));
    double Pv_exact[9*N] __attribute__ ((aligned(32)));
    fill(xv_exact, 2*N, 0.0);
    Particle particles_exact[N];
    init_particles_contigous(particles_exact, xv_exact, Pv_exact, N);
    
    int iwp_exact = 1;
    double G_exact = 0.0;
    Vector3d xtrue_exact;
    data_loader(wp, N_waypoints, V, *Q, dt, N, xtrue_exact, &iwp_exact, &G_exact, particles_exact);

    predict_update_base(wp, N_waypoints, V, *Q, dt, N, xtrue_exact, &iwp_exact, &G_exact, particles_exact);
    

    expect(that % fabs(G - G_exact) < 1.0e-10) << "G";
    expect(that % (iwp - iwp_exact) == 0 ) << "iwp";

    for (int i=0; i<N; i++) {
        expect(that % fabs(particles[i].xv[0] - particles_exact[i].xv[0]) < 1.0e-6) <<particles[i].xv[0] << "x"<< i <<particles_exact[i].xv[0];
        expect(that % fabs(particles[i].xv[1] - particles_exact[i].xv[1]) < 1.0e-6) <<particles[i].xv[1] << "y"<< i <<particles_exact[i].xv[1];
        expect(that % fabs(particles[i].xv[2] - particles_exact[i].xv[2]) < 1.0e-6);
    }

    expect(that % fabs(xtrue[0] - xtrue_exact[0]) < 1.0e-10);
    expect(that % fabs(xtrue[1] - xtrue_exact[1]) < 1.0e-10);
    expect(that % fabs(xtrue[2] - xtrue_exact[2]) < 1.0e-10);

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&predict_update_base)> bench("predict_update Benchmark");

    bench.data_loader = data_loader;
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&predict_update_base, "base", 0.0);
    bench.add_function(&predict_update_old, "functions inplace", 0.0);
    bench.add_function(&predict_update_sine, "sine approximation", 0.0);
    bench.add_function(&predict_update_simd, "basic simd", 0.0);
    bench.add_function(&predict_update_fast_scalar_pipi, "simd + scalar pi_to_pi", 0.0);
    bench.add_function(&predict_update_fast_normal_rand, "simd opt. normal rand", 0.0);
    bench.add_function(&predict_update_fast_plain, "simd opt. fast rand", 0.0);
    bench.add_function(&predict_update_fast, "fastest", 0.0);

    //Run the benchmark: give the inputs of your function in the same order as they are defined.
    set_work(bench, wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles) ;
    bench.run_benchmark(wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles);

    G = M_PI;
    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    set_work(bench, wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles) ;
    bench.run_benchmark(wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles);

    minD = 5.0;
    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    set_work(bench, wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles) ;
    bench.run_benchmark(wp, N_waypoints, V, *Q, dt, N, xtrue, &iwp, &G,particles);

    //bench.destructor_output = false;
    bench.details();

    Benchmark<decltype(&predict_update_base)> bench_scale("predict_update with Particles");
    bench_scale.add_function(&predict_update_base, "base", 0.0);
    bench_scale.add_function(&predict_update_old, "functions inplace", 0.0);
    bench_scale.add_function(&predict_update_sine, "sine approximation", 0.0);
    bench_scale.add_function(&predict_update_simd, "basic simd", 0.0);
    bench_scale.add_function(&predict_update_fast_scalar_pipi, "simd + scalar pi_to_pi", 0.0);
    bench_scale.add_function(&predict_update_fast_normal_rand, "simd opt. normal rand", 0.0);
    bench_scale.add_function(&predict_update_fast_plain, "simd opt. fast rand", 0.0);
    bench_scale.add_function(&predict_update_fast, "fastest", 0.0);
    bench_scale.add_function(&predict_update_fast, "fast", 0.0);

    bench_scale.data_loader = data_loader;
    bench_scale.csv_path = "predict_update_scale_particles.csv";
    bench_scale.csv_output = false;

    int Np = 100;
    for (int i = 0; i< 9; i++) {
        const int Npi = std::pow(2,i) * Np;
        Particle* ps = (Particle*) aligned_alloc(32, Npi * sizeof(Particle));
        double xvi[3*Npi] __attribute__ ((aligned(32)));
        double Pvi[9*Npi] __attribute__ ((aligned(32)));
        fill(xvi, 2*Npi, 0.0);
        init_particles_contigous(ps, xvi, Pvi, Npi);
        set_work(bench_scale, wp, N_waypoints, V, *Q, dt, Npi, xtrue, &iwp, &G,ps);
        bench_scale.run_name = std::to_string(Npi);
        bench_scale.run_benchmark(wp, N_waypoints, V, *Q, dt, Npi, xtrue, &iwp, &G,ps);
        cleanup_members(ps, Npi);
        free(ps);
    }
    bench_scale.details();
    bench_scale.write_csv_details();
    return 0;
}
