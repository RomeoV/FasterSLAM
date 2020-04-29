#include "rdtsc_benchmark.h"
#include <iostream>
#include "add_feature.h"
#include "linalg.h"

int main() {
    // Initialize Input
    Particle p;
    double xv_initial[3] =  {0,0,0};
    initParticle(&p, 10, xv_initial);
    Vector3d pos = {1., 1., acos(4. / 5.)};
    copy(pos, 3, p.xv);
    p.Nfa = 3;

    Vector2d landmarks[2] = {
        {5, 0}, {2, M_PI / 2 - acos(4. / 5.)}};  // local robot frame

    Matrix2d R = {pow(0.1, 2), 0,
                0, pow(1.0 * M_PI / 180, 2)};  // this is from the yglee config

    const size_t N = 1000;

    double(*lambda)(Particle p, Vector2d* landmarks, Matrix2d R) = [](Particle p, Vector2d* landmarks, Matrix2d R){
        double sum=0;
        for (int i = 0; i<N; i++) {
            add_feature(&p, landmarks, 2, R);
        }
        return sum;
    };

    double(*lambda_base)(Particle p, Vector2d* landmarks, Matrix2d R) = [](Particle p, Vector2d* landmarks, Matrix2d R){
        double sum=0;
        for (int i = 0; i<N;i++) {
            add_feature(&p, landmarks, 2, R);
        }
        return sum;
    };

    /*double(*lambda_fmod)(double* angles) = [](double* angles){
        double sum = 0;
        for (int i = 0; i<N;i++) {
            sum+=add_feature(angles[i]);
        }
        return sum;
    };*/

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&add_feature)> bench("add_feature Benchmark");

    double work = 50000.0; // best-case in flops

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&add_feature, "add_feature", work);
    bench.add_function(&add_feature_base, "add_feature_base", work);
    //bench.add_function(&add_feature_fmod, "add_feature_fmod", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(&p, landmarks, 2, R);

    
    Benchmark<decltype(lambda)> bench_lambda("add_feature on array Benchmark");

    // Add lambda functions to aggregate over range of inputs.
    bench_lambda.add_function(lambda, "add_feature", work*N);
    bench_lambda.add_function(lambda_base, "add_feature_base", work*N);
    //bench_lambda.add_function(lambda_fmod, "add_feature_fmod", work*N);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench_lambda.run_benchmark(&p, landmarks, 2, R);


    // Free memory
    // destroy(A);
    // destroy(x);
    // destroy(y);
    delParticleMembers(&p);

    return 0;
}