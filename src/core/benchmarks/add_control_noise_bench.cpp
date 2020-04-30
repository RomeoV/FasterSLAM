#include "rdtsc_benchmark.h"
#include <iostream>
#include "add_control_noise.h"

int main() {
    // Initialize Input
    double V = 10;
    double G = 0.5;
    double Q[4] = {0.1, 0.1, 0.3, 0.5};
    int addnoise = 1;
    double VnGn[2];

    const size_t N = 1000;

    double(*lambda)(double V, double G, double *Q, int addnoise, double *VnGn) = [](double V, double G, double *Q, int addnoise, double *VnGn){
        double sum=0;
        for (int i = 0; i<N; i++) {
            add_control_noise(V,G,Q,addnoise,VnGn);
        }
        return sum;
    };

    double(*lambda_base)(double V, double G, double *Q, int addnoise, double *VnGn) = [](double V, double G, double *Q, int addnoise, double *VnGn){
        double sum=0;
        for (int i = 0; i<N;i++) {
            add_control_noise(V,G,Q,addnoise,VnGn);
        }
        return sum;
    };

    /*double(*lambda_fmod)(double* angles) = [](double* angles){
        double sum = 0;
        for (int i = 0; i<N;i++) {
            sum+=add_control_noise_fmod(angles[i]);
        }
        return sum;
    };*/

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&add_control_noise)> bench("add_control_noise Benchmark");

    double work = 50000.0; // best-case in flops

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&add_control_noise, "add_control_noise", work);
    bench.add_function(&add_control_noise_base, "add_control_noise_base", work);
    //bench.add_function(&add_control_noise_fmod, "pi_to_pi_fmod", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(V,G,Q,addnoise,VnGn);

    
    Benchmark<decltype(lambda)> bench_lambda("add_control_noise on array Benchmark");

    // Add lambda functions to aggregate over range of inputs.
    bench_lambda.add_function(lambda, "add_control_noise", work*N);
    bench_lambda.add_function(lambda_base, "add_control_noise_base", work*N);
    //bench_lambda.add_function(lambda_fmod, "add_control_noise_fmod", work*N);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench_lambda.run_benchmark(V,G,Q,addnoise,VnGn);


    // Free memory
    // destroy(A);
    // destroy(x);
    // destroy(y);

    return 0;
}
