#include "rdtsc_benchmark.h"
#include <iostream>
#include "add_feature.h"
#include "linalg.h"

void data_loader(Particle *p, Vector2d *z, size_t N_z, double *R) {
    // Initialize Input
    double xv_initial[3] =  {0,0,0};
    initParticle(p, 300000, xv_initial);
    Vector3d pos = {1., 1., acos(4. / 5.)};
    copy(pos, 3, p->xv);
    p->Nfa = 3;
}

int main() {
    
    Particle p;
    Vector2d landmarks[2] = {
        {5, 0}, {2, M_PI / 2 - acos(4. / 5.)}};  // local robot frame

    Matrix2d R = {pow(0.1, 2), 0,
                0, pow(1.0 * M_PI / 180, 2)};  // this is from the yglee config

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&add_feature)> bench("add_feature Benchmark");

    // plus function calls, 2x matrix mul
    double work = 8.0; // best-case in flops

    bench.data_loader = data_loader;
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&add_feature, "add_feature", work);
    bench.add_function(&add_feature_base, "add_feature_base", work);
    //bench.add_function(&add_feature_fmod, "add_feature_fmod", work);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(&p, landmarks, 2, R);

    // Free memory
    delParticleMembers(&p);

    return 0;
}