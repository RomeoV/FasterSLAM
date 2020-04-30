#include "rdtsc_benchmark.h"
#include <iostream>
#include "feature_update.h"

int main() {
    // Initialize Input
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
    size_t N_idf = 3;

    Vector2d z[2] = {{25.77444969441645, -1.4733774573801159}, {25.276107769232738, 0.13836675004531551}};
    const int zlen = 2;
    double R[4] = {0.010000000000000002, 0, 0, 0.00030461741909055634};
    const int addnoise = 1;

    const size_t N = 1000;

    double(*lambda)(Particle* particle, Vector2d* z, int* idf, size_t N_idf, Matrix2d R) = [](Particle* particle, Vector2d* z, int* idf, size_t N_idf, Matrix2d R){
        double sum=0;
        for (int i = 0; i<N; i++) {
            feature_update(particle, z, idf, N_idf, R);
        }
        return sum;
    };

    double(*lambda_base)(Particle* particle, Vector2d* z, int* idf, size_t N_idf, Matrix2d R) = [](Particle* particle, Vector2d* z, int* idf, size_t N_idf, Matrix2d R){
        double sum=0;
        for (int i = 0; i<N;i++) {
            feature_update(particle, z, idf, N_idf, R);
        }
        return sum;
    };

    /*double(*lambda_fmod)(double* angles) = [](double* angles){
        double sum = 0;
        for (int i = 0; i<N;i++) {
            sum+=feature_update(angles[i]);
        }
        return sum;
    };*/

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&feature_update)> bench("feature_update Benchmark");

    double work = 50000.0; // best-case in flops

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&feature_update, "feature_update", work);
    bench.add_function(&feature_update_base, "feature_update_base", work);
    //bench.add_function(&feature_update_fmod, "feature_update_fmod", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(particle, z, idf, N_idf, R);

    
    Benchmark<decltype(lambda)> bench_lambda("feature_update on array Benchmark");

    // Add lambda functions to aggregate over range of inputs.
    bench_lambda.add_function(lambda, "feature_update", work*N);
    bench_lambda.add_function(lambda_base, "afeature_update_base", work*N);
    //bench_lambda.add_function(lambda_fmod, "feature_update_fmod", work*N);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench_lambda.run_benchmark(particle, z, idf, N_idf, R);


    // Free memory
    // destroy(A);
    // destroy(x);
    // destroy(y);

    return 0;
}
