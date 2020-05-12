#include "rdtsc_benchmark.h"
#include <iostream>
#include "particle.h"

int main() {
    // // Initialize Input
    // double *A, *x, *y;
    // const size_t N = 1000;

    // double(*lambda)(double* angles) = [](double* angles){
    //     double sum=0;
    //     for (int i = 0; i<N; i++) {
    //         //sum+=copyParticle(angles[i]);
    //     }
    //     return sum;
    // };

    // double(*lambda_base)(double* angles) = [](double* angles){
    //     double sum=0;
    //     for (int i = 0; i<N;i++) {
    //         //sum+=copyParticle_base(angles[i]);
    //     }
    //     return sum;
    // };

    // /*double(*lambda_fmod)(double* angles) = [](double* angles){
    //     double sum = 0;
    //     for (int i = 0; i<N;i++) {
    //         sum+=copyParticle_fmod(angles[i]);
    //     }
    //     return sum;
    // };*/

    // // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    // //Benchmark<decltype(&copyParticle)> bench("copyParticle Benchmark");

    // double work = 50000.0; // best-case in flops

    // // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // // First function should always be the base case you want to benchmark against!
    // //bench.add_function(&copyParticle, "copyParticle", work);
    // //bench.add_function(&copyParticle_base, "copyParticle_base", work);
    // //bench.add_function(&copyParticle_fmod, "pi_to_pi_fmod", work);

    // //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    // //bench.run_benchmark(angles[0]);

    
    // Benchmark<decltype(lambda)> bench_lambda("copyParticle on array Benchmark");

    // Add lambda functions to aggregate over range of inputs.
    //bench_lambda.add_function(lambda, "copyParticle", work*N);
    //bench_lambda.add_function(lambda_base, "copyParticle_base", work*N);
    //bench_lambda.add_function(lambda_fmod, "copyParticle_fmod", work*N);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    //bench_lambda.run_benchmark(angles);


    // Free memory
    // destroy(A);
    // destroy(x);
    // destroy(y);

    return 0;
}
