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

#include "observe_update.h"
#include "predict_update.h"

#include "get_observations.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

# define Nfa_start 0

void data_loader(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    srand(0);
    double minD = 0.1;
    double rateG = 1.0;
    double maxG = 1.5;
    double dt = 1.0;
    double V = 3.0;
    setup_initial_Q_R();
    double xv_initial[3] = {0.0,0.0,0.0};
    int N_waypoints = 3;
    double wp[6] = {0,0,1,1,2,2};
    int iwp = 1;
    double G = 0.0;
    *Nf_visible = 0;
    //double weights[NPARTICLES] __attribute__ ((aligned(32)));
    setup_landmarks(&ftag, &da_table, N_features);
    for (int i = 0; i<Nfa_start; i++) {
        double _xf[2];
        double _Pf[4];
        // Add here start values to each particle to not test the add_feature part!
    }
    for (int i = 0; i<NPARTICLES; i++) {
        particles[i].Nfa = 0;
        weights[i] = 1.0/NPARTICLES;
        particles[i].xv[0] = xv_initial[0];
        particles[i].xv[1] = xv_initial[1];
        particles[i].xv[2] = xv_initial[2];
    }

    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

            //z is the range and bearing of the observed landmark
    
    get_observations_base(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    
    predict_update_base(wp, N_waypoints, V, *Q, dt, NPARTICLES, xv_initial, &iwp, &G,particles);
    observe_update_base(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf, zn, particles, weights);

    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

            //z is the range and bearing of the observed landmark
    
    get_observations_base(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
}

void init_particles_contigous(Particle* particle, double* _xv, double* _Pv, double* weight, const size_t N, const size_t Nf) {
    for (int i = 0; i<N; i++) {
        particle[i].xv = _xv+3*i;
        particle[i].Pv = _Pv+9*i;
        weight[i] = 1.0/N;
        particle[i].w = weight+i;
        initParticle_prealloc(particle+i, Nf, _xv+3*i);

        double xf[2];
        double Pf[4];
    } 
}

// I will try to add this as smooth as possible to the benchmark, but for now do this
void set_work(Benchmark<decltype(&observe_update)>& bench, 
                double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    data_loader(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf, zn, particles, weights);
    bench.funcFlops[0] = observe_update_base_flops(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf, zn, particles, weights);
    data_loader(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf, zn, particles, weights);
    bench.funcBytes[0] = 8 * observe_update_base_memory(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf, zn, particles, weights);
    data_loader(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf, zn, particles, weights);

    double active_flops = observe_update_flops(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf, zn, particles, weights);
        data_loader(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf, zn, particles, weights);
    double active_memory = 8* observe_update_memory(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf, zn, particles, weights);
        data_loader(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf, zn, particles, weights);

    for (int i = 1; i < bench.numFuncs; i++) {
        bench.funcFlops[i] = active_flops;
        bench.funcBytes[i] = active_memory;
    }
}


void setup(Particle* particles, double* weights, const int Nf, Vector3d xtrue, 
                    int** ftag, int** da_table, Vector2d** z, Vector2d** zf, 
                    Vector2d** zn, int** idf, int** ftag_visible, double* xv, double* Pv) {
    init_particles_contigous(particles, xv, Pv, weights, NPARTICLES, Nf);

    // setup_initial_Q_R();  // modifies global variables
    setup_landmarks(ftag, da_table, Nf);
    setup_measurements(z, zf, zn, idf, ftag_visible, Nf);
}

void cleanup_members(Particle* particles, int N) {
    for(int i = 0; i<N; i++) {
        delParticleMembers_prealloc(particles+i);
    }
}

int main() {
    SWITCH_PREDICT_NOISE=1;
    SWITCH_CONTROL_NOISE=1;
    const int Nf = 20;
    double lm[2*Nf];
    const int N = 100; // Cannot set to NPARTICLES...
    
    setup_initial_Q_R();  // Remove this soon!!!
    double xtrue[3] = {0,0,0};
    int N_waypoints = 3;
    double wp[6] = {0,0,1,1,2,2};
    
    
    srand(0);
    fill_rand(lm, 2*Nf, -MAX_RANGE,MAX_RANGE);

    srand(0);
    size_t Nf_visible;
    Vector2d *z;
    Vector2d *zf;
    Vector2d *zn;
    int *idf, *ftag_visible, *ftag, *da_table;
    double xvs[3*N] __attribute__ ((aligned(32)));
    double Pvs[9*N] __attribute__ ((aligned(32)));
    fill(xvs, 3*N, 0.0);
    Particle particles[N]  __attribute__ ((aligned(32)));
    double weights[N]  __attribute__ ((aligned(32)));
    

    setup(particles, weights, Nf, xtrue, &ftag, &da_table, &z, &zf, &zn, &idf, &ftag_visible, xvs, Pvs);

    for (size_t i = 0; i < Nf; i++) {
        ftag_visible[i] = ftag[i];
    }
    get_observations_base(xtrue, MAX_RANGE, lm, Nf, ftag_visible, &Nf_visible, z);
    
    observe_update_fast(lm,Nf, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);

    observe_update_fast(lm,Nf, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);

    
    srand(0);
    size_t Nf_visible_exact;
    Vector2d *z_exact; 
    Vector2d *zf_exact;
    Vector2d *zn_exact;
    int *idf_exact, *ftag_visible_exact, *ftag_exact, *da_table_exact;
    double xv_exact[3*N] __attribute__ ((aligned(32)));
    double Pv_exact[9*N] __attribute__ ((aligned(32)));
    fill(xv_exact, 3*N, 0.0);
    Particle particles_exact[N];
    double weights_exact[N];
    
    setup(particles_exact, weights_exact, Nf, xtrue, &ftag_exact, &da_table_exact, 
        &z_exact, &zf_exact, &zn_exact, &idf_exact, &ftag_visible_exact, xv_exact, Pv_exact);

    for (size_t i = 0; i < Nf; i++) {
        ftag_visible_exact[i] = ftag_exact[i];
    }
    get_observations_base(xtrue, MAX_RANGE, lm, Nf, ftag_visible_exact, &Nf_visible_exact, z_exact);
    observe_update_base(lm,Nf, xtrue, *R, ftag_exact, 
            da_table_exact, ftag_visible_exact, z_exact, &Nf_visible_exact, zf_exact, idf_exact, 
            zn_exact, particles_exact, weights_exact);

    observe_update_base(lm,Nf, xtrue, *R, ftag_exact, 
            da_table_exact, ftag_visible_exact, z_exact, &Nf_visible_exact, zf_exact, idf_exact, 
            zn_exact, particles_exact, weights_exact);

    //TEST
    double eps = 1e-10;
    auto is_close = [&](double lhs, double rhs) -> bool {return std::abs(lhs - rhs) < eps;};
    expect(that % Nf_visible == Nf_visible_exact) << "Nf_visible";
    // print(weights, 100,1);
    // print(weights_exact, 100,1);

    for (int i =0; i<N; i++) {
        expect(that % fabs(weights[i]- weights_exact[i]) < eps)<< "weights";
        for (int j = 0; j<3; j++) {
            expect(that % fabs(particles[i].xv[j] - particles_exact[i].xv[j]) < eps) << "xv" << j;
        }

        for (int j =0; j<9; j++) {
            expect(that % fabs(particles[i].Pv[j] - particles_exact[i].Pv[j]) < eps) << "Pv" << j;
        }
        assert(particles[i].Nfa == particles_exact[i].Nfa);
        for (int j =0; j; j++) {
            expect(that % fabs(particles[i].xf[2*j] - particles_exact[i].xf[2*j]) < eps) << "xf0" << j;
            expect(that % fabs(particles[i].xf[2*j+1] - particles_exact[i].xf[2*j+1]) < eps) << "xf1" << j;
            expect(that % fabs(particles[i].Pf[4*j+0] - particles_exact[i].Pf[4*j+0]) < eps) << "Pf0" << j;
            expect(that % fabs(particles[i].Pf[4*j+1] - particles_exact[i].Pf[4*j+1]) < eps) << "Pf1" << j;
            expect(that % fabs(particles[i].Pf[4*j+2] - particles_exact[i].Pf[4*j+2]) < eps) << "Pf2" << j;
            expect(that % fabs(particles[i].Pf[4*j+3] - particles_exact[i].Pf[4*j+3]) < eps) << "Pf3" << j;
        }
    }
    for (int i = 0; i<Nf; i++) {
        for (int j = 0; j<2; j++) {
            expect(that % fabs(z[i][j] - z_exact[i][j] ) < eps) <<i << "z" << j;
            expect(that % fabs(zf[i][j]  - zf_exact[i][j] ) < eps) << "zf" << j;
            expect(that % fabs(zn[i][j]  - zn_exact[i][j] ) < eps) << "zn" << j;
        }
        expect(that % idf[i] == idf_exact[i])  << "idf" << i;
        expect(that % ftag_visible[i] == ftag_visible_exact[i])  << "ftag_visible" << i;
        expect(that % ftag[i] == ftag_exact[i])  << "ftag" << i;
        expect(that % da_table[i] == da_table_exact[i])  << "da_table" << i;
    }

    // auto ymm0 = _mm256_set_pd(3,2,1,0);
    // auto ymm1 = _mm256_set_pd(7,6,5,4);

    // print256d(_mm256_permute2f128_pd(ymm0,ymm1,0b0001));
    // print256d(_mm256_permute2f128_pd(ymm0,ymm1,0b00100000));


    // print256d(exp_avx2_pd(_mm256_set_pd(3.0,2.0,-63.4411f,-207.536f)));

    // std::cout<<exp(-207.536)<<", "<<exp(-63.4411)<<", "<<exp(2.0)<<", "<<exp(3.0)<<std::endl;

    Benchmark<decltype(&observe_update_base)> bench("observe_update Benchmark");

    bench.data_loader=data_loader;

    bench.controls.NUM_RUNS = 3;

    bench.add_function(&observe_update_base, "observe_update_base", 0.0);
    bench.add_function(&observe_update_active, "observe_update_active", 0.0);

    //Top-Level-functions inlined
    bench.add_function(&observe_update_inplace, "observe_update_inplace", 0.0);

    // Reduced ops
    bench.add_function(&observe_update_simplified, "observe_update_simplified", 0.0);
    bench.add_function(&observe_update_unrolled4x, "observe_update_unrolled4x", 0.0);
    bench.add_function(&observe_update_active, "observe_update_unrolled4x_plain", 0.0);

    bench.add_function(&observe_update_simpleavx, "observe_update_simpleavx", 0.0);
    bench.add_function(&observe_update_fast_KF_comp_not_unrolled, "observe_update_fast_KF_comp_not_unrolled", 0.0);
    bench.add_function(&observe_update_fast_partial_avx_slower_vTMv, "observe_update_fast_improved_vTMv", 0.0);
    bench.add_function(&observe_update_fast, "observe_update_fast_jacobians", 0.0);
#ifndef KF_YGLEE
    bench.add_function(&observe_update_fast_fullavx, "observe_update_all_unrolled", 0.0);
#endif
    bench.add_function(&observe_update, "observe_update", 0.0);

    set_work(bench, lm,Nf, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn_exact, particles, weights);
    bench.run_benchmark(lm,Nf, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn_exact, particles, weights);

    Benchmark<decltype(&observe_update_base)> bench_scale("observe_update Scale Particles");

    bench_scale.data_loader=data_loader;
    bench_scale.controls.NUM_RUNS = 3;
    bench_scale.add_function(&observe_update_base, "observe_update_base", 0.0);
    bench_scale.add_function(&observe_update_simplified, "observe_update_simplified", 0.0);
    bench_scale.add_function(&observe_update_unrolled4x, "observe_update_unrolled4x", 0.0);
    bench_scale.add_function(&observe_update_active, "observe_update_unrolled4x_plain", 0.0);

    bench_scale.add_function(&observe_update_simpleavx, "observe_update_simpleavx", 0.0);
    bench_scale.add_function(&observe_update_fast_KF_comp_not_unrolled, "observe_update_fast_KF_comp_not_unrolled", 0.0);
    bench_scale.add_function(&observe_update_fast_partial_avx_slower_vTMv, "observe_update_fast_improved_vTMv", 0.0);
    bench_scale.add_function(&observe_update_fast, "observe_update_fast_jacobians", 0.0);
#ifndef KF_YGLEE
    bench_scale.add_function(&observe_update_fast_fullavx, "observe_update_fast", 0.0);
#endif    

    bench_scale.data_loader = data_loader;
    bench_scale.csv_path = "observe_update_scale_particles.csv";
    bench_scale.csv_output = false;

    int Np = 100;
    for (int i = 0; i< 9; i++) {
        NPARTICLES = std::pow(2,i) * Np;
        std::cout<< "Benchmarking N="<<NPARTICLES<<" Particles..."<<std::endl;
         double* ws = (double*) aligned_alloc(32, NPARTICLES * sizeof(double));
        Particle* ps = (Particle*) aligned_alloc(32, NPARTICLES* sizeof(Particle));
        double xvi[3*NPARTICLES] __attribute__ ((aligned(32)));
        double Pvi[9*NPARTICLES] __attribute__ ((aligned(32)));
        fill(xvi, 2*NPARTICLES, 0.0);
        init_particles_contigous(ps, xvi, Pvi, ws, NPARTICLES, Nf);
        set_work(bench_scale, lm,Nf, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn_exact, ps, ws);
        bench_scale.run_name = std::to_string(NPARTICLES);
        bench_scale.run_benchmark(lm,Nf, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn_exact, ps, ws);

        cleanup_members(ps, NPARTICLES);
        free(ps);
        free(ws);
    }

    bench_scale.details();
    bench_scale.write_csv_details();
    return 0;
}