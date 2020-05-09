#include "rdtsc_benchmark.h"
#include <iostream>
#include "add_feature.h"
#include "linalg.h"

Vector2d landmarks[2] = {
    {5, 0}, {2, M_PI / 2 - acos(4. / 5.)}};  // local robot frame

Matrix2d R = {pow(0.1, 2), 0,
            0, pow(1.0 * M_PI / 180, 2)};  // this is from the yglee config

void add_features_bench(Particle* p, Vector2d* landmarks, Matrix2d R) {
    for (size_t i = 0; i < 100; i++) {
        p->Nfa = 0;
        add_feature(p, landmarks, 2, R);
    }
}

void add_features_base_bench(Particle* p, Vector2d* landmarks, Matrix2d R) {
    for (size_t i = 0; i < 100; i++) {
        p->Nfa = 0;
        add_feature_base(p, landmarks, 2, R);
    }
}

int main() {
    // Initialize Input
    Particle p;
    double xv_initial[3] =  {0,0,0};
    initParticle(&p, 3, xv_initial);
    Vector3d pos = {1., 1., acos(4. / 5.)};
    copy(pos, 3, p.xv);
    p.Nfa = 0;

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&add_features_bench)> bench("add_feature Benchmark");

    double work = 50000.0; // best-case in flops

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&add_features_bench, "add_feature", work);
    bench.add_function(&add_features_base_bench, "add_feature_base", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(&p, landmarks, R);

    
    delParticleMembers(&p);

    return 0;
}