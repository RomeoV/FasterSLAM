#include "rdtsc_benchmark.h"
#include <iostream>
#include "compute_jacobians.h"


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
    int idf[3] = {0,0,0};
    int N_z = 3;
    Matrix2d R = {1,0,0,1};

    // outputs
    Vector2d zp[2*3] = {0,0, 0,0, 0,0}; // measurement (range, bearing)
    Matrix23d Hv[6*3] = {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0}; // jacobians of function h (deriv of h wrt pose)
    Matrix2d Hf[4*3] = {0,0,0,0, 0,0,0,0, 0,0,0,0}; // jacobians of function h (deriv of h wrt mean)
    Matrix2d Sf[4*3] = {0,0,0,0, 0,0,0,0, 0,0,0,0}; // Measurement covariance of feature observation given the vehicle.

    const size_t N = 1000;

    double(*lambda)(Particle* particle, int* idf, int N_z, Matrix2d R, Vector2d* zp, Matrix23d* Hv, Matrix2d* Hf, Matrix2d* Sf) = [](Particle* particle, int* idf, int N_z, Matrix2d R, Vector2d* zp, Matrix23d* Hv, Matrix2d* Hf, Matrix2d* Sf){
        double sum=0;
        for (int i = 0; i<N; i++) {
            compute_jacobians(particle, idf, N_z, R, zp, Hv, Hf, Sf);
        }
        return sum;
    };

    double(*lambda_base)(Particle* particle, int* idf, int N_z, Matrix2d R, Vector2d* zp, Matrix23d* Hv, Matrix2d* Hf, Matrix2d* Sf) = [](Particle* particle, int* idf, int N_z, Matrix2d R, Vector2d* zp, Matrix23d* Hv, Matrix2d* Hf, Matrix2d* Sf){
        double sum=0;
        for (int i = 0; i<N;i++) {
            compute_jacobians(particle, idf, N_z, R, zp, Hv, Hf, Sf);
        }
        return sum;
    };

    /*double(*lambda_fmod)(double* angles) = [](double* angles){
        double sum = 0;
        for (int i = 0; i<N;i++) {
            sum+=compute_jacobians(angles[i]);
        }
        return sum;
    };*/

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&compute_jacobians)> bench("compute_jacobians Benchmark");

    double work = 50000.0; // best-case in flops

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&compute_jacobians_base, "compute_jacobians_base", work);
    bench.add_function(&compute_jacobians, "compute_jacobians", work);
    //bench.add_function(&compute_jacobians_fmod, "compute_jacobians_fmod", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(particle, idf, N_z, R, zp, Hv, Hf, Sf);

    
    /*Benchmark<decltype(lambda)> bench_lambda("compute_jacobians on array Benchmark");

    // Add lambda functions to aggregate over range of inputs.
    bench_lambda.add_function(lambda_base, "compute_jacobians_base", work*N);
    bench_lambda.add_function(lambda, "compute_jacobians", work*N);
    //bench_lambda.add_function(lambda_fmod, "compute_jacobians_fmod", work*N);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench_lambda.run_benchmark(particle, idf, N_z, R, zp, Hv, Hf, Sf);*/


    // Free memory
    // destroy(A);
    // destroy(x);
    // destroy(y);
    //! Delete Particle
    delParticleMembersAndFreePtr(particle);

    return 0;
}
