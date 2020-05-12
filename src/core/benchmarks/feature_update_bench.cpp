#include "rdtsc_benchmark.h"
#include <iostream>
#include "feature_update.h"

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
    int idf[3] = {0,0,0};
    size_t N_idf = 3;

    Vector2d z[2] = {{25.77444969441645, -1.4733774573801159}, {25.276107769232738, 0.13836675004531551}};
    const int zlen = 2;
    double R[4] = {0.010000000000000002, 0, 0, 0.00030461741909055634};
    const int addnoise = 1;

    const size_t N = 1000;

    /*void(*lambda)(Particle *particle, Vector2d *z, int *idf, size_t N_idf, double *R) = [](Particle *particle, Vector2d *z, int *idf, size_t N_idf, double *R){
        for (int i = 0; i<N; i++) {
            feature_update(particle, z, idf, N_idf, R);
        }
    };

    void(*lambda_base)(Particle *particle, Vector2d *z, int *idf, size_t N_idf, double *R) = [](Particle *particle, Vector2d *z, int *idf, size_t N_idf, double *R){
        for (int i = 0; i<N;i++) {
            feature_update_base(particle, z, idf, N_idf, R);
        }
    };*/

    /*double(*lambda_fmod)(double* angles) = [](double* angles){
        double sum = 0;
        for (int i = 0; i<N;i++) {
            sum+=feature_update_fmod(angles[i]);
        }
        return sum;
    };*/

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&feature_update)> bench("feature_update Benchmark");

    double work = 500.0; // best-case in flops

    bench.data_loader = data_loader;
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&feature_update, "feature_update", work);
    bench.add_function(&feature_update_base, "feature_update_base", work);
    //bench.add_function(&feature_update_fmod, "feature_update_fmod", work);

    Vector2d zp[3];
    Matrix23d Hv[3];
    Matrix2d Hf[3];
    Matrix2d Sf[3];
    compute_jacobians_base(&p, idf, 3, R, zp, Hv, Hf, Sf);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(&p, z, idf, N_idf, R, zp, Hv, Hf, Sf);

    /*Benchmark<decltype(lambda)> bench_lambda("feature_update on array Benchmark");

    bench_lambda.data_loader = data_loader;
    // Add lambda functions to aggregate over range of inputs.
    bench_lambda.add_function(lambda, "feature_update", work*N);
    bench_lambda.add_function(lambda_base, "afeature_update_base", work*N);
    //bench_lambda.add_function(lambda_fmod, "feature_update_fmod", work*N);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench_lambda.run_benchmark(&p, z, idf, N_idf, R);*/


    // Free memory
    delParticleMembers(&p);

    return 0;
}
