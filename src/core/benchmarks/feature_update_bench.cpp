#include "rdtsc_benchmark.h"
#include <iostream>
#include "feature_update.h"
#include "linalg.h"

int main() {
    
    Particle p;
    int idf[3] = {0,0,0};
    size_t N_idf = 3;
    double xv_initial[3] =  {0,0,0};
    initParticle(&p, 300000, xv_initial);
    Vector3d pos = {1., 1., acos(4. / 5.)};
    copy(pos, 3, p.xv);
    p.Nfa = 3;
    Vector2d xf[3] = {{1,0.1},{1,0.2},{1,0.3}}; 
    for(int i=0; i<3; i++){ //! 2d means of EKF in cartesian world coordinates
        set_xfi(&p, xf[i], i);
    }

    Matrix2d Pf[3] = {{1,0,1,0},{1,0,1,0},{1,0,1,0}}; //! covariance matrices for EKF in polar robot coordinates
        for(int i=0; i<3; i++){ //! covariance matrices for EKF (=Extended Kalman Filter) in polar robot coordinates
        set_Pfi(&p, Pf[i], i);
    }

    Vector2d z[2] = {{25.77444969441645, -1.4733774573801159}, {25.276107769232738, 0.13836675004531551}};
    const int zlen = 2;
    double R[4] = {0.010000000000000002, 0, 0, 0.00030461741909055634};
    const int addnoise = 1;

    const size_t N = 1000;

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&feature_update)> bench("feature_update Benchmark");

    double work = 500.0; // best-case in flops

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&feature_update, "feature_update", work);
    bench.add_function(&feature_update_base, "feature_update_base", work);

    Vector2d zp[3];
    Matrix23d Hv[3];
    Matrix2d Hf[3];
    Matrix2d Sf[3];
    compute_jacobians_base(&p, idf, 3, R, zp, Hv, Hf, Sf);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(&p, z, idf, N_idf, R, zp, Hv, Hf, Sf);

    // Free memory
    delParticleMembers(&p);

    return 0;
}
