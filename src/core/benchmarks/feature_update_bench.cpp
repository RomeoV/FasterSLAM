#include "rdtsc_benchmark.h"
#include <iostream>
#include "feature_update.h"
#include "particle.h"

#include "ut.hpp"
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

void data_loader(Particle *p, Vector2d *z, int *idf, size_t N_idf, double *R, Vector2d *zp, Matrix23d *Hv, Matrix2d *Hf, Matrix2d *Sf) {
    // Initialize Input
    Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
    initParticle(p, 300000, xv);
    Vector2d xf[3] = {{1,0.1},{1,0.2},{1,0.3}}; 
    for(int i=0; i<3; i++){ //! 2d means of EKF in cartesian world coordinates
        set_xfi(p, xf[i], i);
    }
    Matrix2d Pf[3] = {{1,0,1,0},{1,0,1,0},{1,0,1,0}}; //! covariance matrices for EKF in polar robot coordinates
        for(int i=0; i<3; i++){ //! covariance matrices for EKF (=Extended Kalman Filter) in polar robot coordinates
        set_Pfi(p, Pf[i], i);
    }
}

int main() {
    
    Particle p;
    double xv_initial[3] = {0,0,0};
    initParticle(&p, 0, xv_initial);
    int idf[3] = {0,0,0};
    size_t N_idf = 3;

    Vector2d z[2] = {{25.77444969441645, -1.4733774573801159}, {25.276107769232738, 0.13836675004531551}};
    const int zlen = 2;
    double R[4] = {0.010000000000000002, 0, 0, 0.00030461741909055634};
    const int addnoise = 1;

    Vector2d zp[3] __attribute__((aligned(32)));
    Matrix23d Hv[3] __attribute__((aligned(32)));
    Matrix2d Hf[3] __attribute__((aligned(32)));
    Matrix2d Sf[3] __attribute__((aligned(32)));
    compute_jacobians_base(&p, idf, 3, R, zp, Hv, Hf, Sf);

    // sanity check
    Particle p_base, p_new;
    double xv[3] = {1., 1., acos(4. / 5.)};
    initParticle(&p_base, 3, xv);
    (&p_base)->Nfa = 3;
    initParticle(&p_new, 3, xv);
    (&p_new)->Nfa = 3;
    // need to set before updating
    Vector2d xf[3] = {{1,0.1},{1,0.2},{1,0.3}}; 
    for(int i=0; i<3; i++){ //! 2d means of EKF in cartesian world coordinates
        set_xfi(&p_base, xf[i], i);
        set_xfi(&p_new, xf[i], i);
    }
    Matrix2d Pf[3] = {{1,0,1,0},{1,0,1,0},{1,0,1,0}}; //! covariance matrices for EKF in polar robot coordinates
    for(int i=0; i<3; i++){ //! covariance matrices for EKF (=Extended Kalman Filter) in polar robot coordinates
        set_Pfi(&p_base, Pf[i], i);
        set_Pfi(&p_new, Pf[i], i);
    }
    
    // not necessary to set xf and Pf 
    feature_update_base(&p_base, z, idf, N_idf, R, zp, Hv, Hf, Sf);
    feature_update(&p_new, z, idf, N_idf, R, zp, Hv, Hf, Sf);

    // check one feature
    for(int i=0; i<2; i++){
        double xfi_new = p_new.xf[2*2+i];
        double xfi_base = p_base.xf[2*2+i];
        expect(that % fabs(xfi_new - xfi_base) < 1.0e-10);
    } 

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&feature_update)> bench("feature_update Benchmark");
    
    double work = 500.0; // best-case in flops

    bench.data_loader = data_loader;
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&feature_update_base, "feature_update_base", work);
    bench.add_function(&feature_update, "feature_update", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(&p, z, idf, N_idf, R, zp, Hv, Hf, Sf);


    // Free memory
    delParticleMembers(&p);
    delParticleMembers(&p_new);
    delParticleMembers(&p_base);

    return 0;
}
