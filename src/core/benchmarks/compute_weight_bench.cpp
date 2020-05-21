#include "rdtsc_benchmark.h"
#include <iostream>
#include "compute_weight.h"
#include "compute_jacobians.h"

#include "ut.hpp"
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

int main() {
    // Initialize Input
    // prepare particle
    Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
    Particle* particle = newParticle(5, xv);
    Vector2d xf[3] = {{1,0.1},{1,0.2},{1,0.3}}; 
    for(int i=0; i<3; i++){ //! 2d means of EKF in cartesian world coordinates
    set_xfi(particle, xf[i], i);
    }
    Matrix2d Pf[3] = {{1,0,1,0},{1,0,1,0},{1,0,1,0}}; //! covariance matrices for EKF in polar robot coordinates
    for(int i=0; i<3; i++){ //! covariance matrices for EKF (=Extended Kalman Filter) in polar robot coordinates
    set_Pfi(particle, Pf[i], i);
    }

    Vector2d z[3] = {{0,0}, {0,0}, {0,0}};  // vector of map features
    size_t N_z = 3; // number of features.
    int idf[3] = {0,0,0};  // vector of map indices
    Matrix2d R = {1,0,0,1};   // matrix of observation noises

    const size_t N = 1000;

    Vector2d zp[3];
    Matrix23d Hv[3];
    Matrix2d Hf[3];
    Matrix2d Sf[3];
    compute_jacobians_base(particle, idf, 3, R, zp, Hv, Hf, Sf);

    // sanity check 
    double weight_base = compute_weight(particle, z, N_z, idf, R, zp, Hv, Hf, Sf);
    double weight_new = compute_weight(particle, z, N_z, idf, R, zp, Hv, Hf, Sf);
    expect(that % fabs(weight_base - weight_new) < 1.0e-10) << "weight";

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&compute_weight)> bench("compute_weight Benchmark");

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&compute_weight_base, "compute_weight_base", 0.0);
    bench.funcFlops[0] = compute_weight_base_flops(particle, z, N_z, idf, R, zp, Hv, Hf, Sf);
    bench.funcBytes[0] = compute_weight_base_memory(particle, z, N_z, idf, R, zp, Hv, Hf, Sf);
    bench.add_function(&compute_weight, "compute_weight", 0.0);
    bench.funcFlops[1] = compute_weight_active_flops(particle, z, N_z, idf, R, zp, Hv, Hf, Sf);
    bench.funcBytes[1] = compute_weight_active_memory(particle, z, N_z, idf, R, zp, Hv, Hf, Sf);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(particle, z, N_z, idf, R, zp, Hv, Hf, Sf);

    //! Delete Particle
    delParticleMembersAndFreePtr(particle);

    return 0;
}
