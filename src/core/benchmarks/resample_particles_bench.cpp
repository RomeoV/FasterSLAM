#include "rdtsc_benchmark.h"
#include <iostream>
#include "resample_particles.h"

int main() {
    // // Initialize Input
    // double weights[3] = {2./3, 1./3, 0.};
    // const size_t Nf = 5;

    // Particle particles[3];

    // Vector3d zeros = {0.,0.,0.};
    // particles[0].w = &weights[0];
    // particles[0].Nf = Nf;

    // Vector3d ones = {1.,1.,1.};
    // particles[1].w = &weights[1];
    // particles[1].Nf = Nf;

    // Vector3d twos = {2.,2.,2.};
    // particles[2].w = &weights[2];
    // particles[2].Nf = Nf;

    // for (size_t i = 0; i < 3; i++) {
    //     auto xv = std::vector{zeros, ones, twos};
    //     initParticle(&particles[i], 5, xv[i]);
    //     particles[i].Nfa = 3;
    //     for (size_t el = 0; el < 2*3; el++) {
    //         particles[i].xf[el] = i;
    //     }
    // }

    // const size_t N = 1000;

    // double(*lambda)(Particle* particles, double* weights) = [](Particle* particles, double* weights){
    //     double sum=0;
    //     for (int i = 0; i<N; i++) {
    //        resample_particles(particles, 3, weights, 2, 1);
    //     }
    //     return sum;
    // };

    // double(*lambda_base)(Particle* particles, double* weights) = [](Particle* particles, double* weights){
    //     double sum=0;
    //     for (int i = 0; i<N;i++) {
    //        resample_particles(particles, 3, weights, 2, 1);
    //     }
    //     return sum;
    // };

    // /*double(*lambda_fmod)(double* angles) = [](double* angles){
    //     double sum = 0;
    //     for (int i = 0; i<N;i++) {
    //         sum+=resample_particles(angles[i]);
    //     }
    //     return sum;
    // };*/

    // // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    // Benchmark<decltype(&resample_particles)> bench("resample_particles Benchmark");

    // double work = 50000.0; // best-case in flops

    // // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // // First function should always be the base case you want to benchmark against!
    // bench.add_function(&resample_particles, "resample_particles", work);
    // bench.add_function(&resample_particles_base, "resample_particles_base", work);
    // //bench.add_function(&resample_particles_fmod, "resample_particles_fmod", work);

    // //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    // bench.run_benchmark(particles, weights);

    
    // Benchmark<decltype(lambda)> bench_lambda("resample_particles on array Benchmark");

    // // Add lambda functions to aggregate over range of inputs.
    // bench_lambda.add_function(lambda, "resample_particles", work*N);
    // bench_lambda.add_function(lambda_base, "resample_particles_base", work*N);
    // //bench_lambda.add_function(lambda_fmod, "resample_particles_fmod", work*N);

    // //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    // bench_lambda.run_benchmark(particles, weights);


    // // Free memory
    // // destroy(A);
    // // destroy(x);
    // // destroy(y);

    return 0;
}