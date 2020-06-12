#include "rdtsc_benchmark.h"
#include <iostream>
#include "resample_particles.h"
#include "linalg.h"


void data_loader(Particle* particles, size_t N, double* weights,int Nmin, int doresample) {
    srand(0);
    fill_rand(weights, N, 0.0001,1.0); //We normalize anyway in resample
}


void setup_particles(Particle* particles, double* weights, const int N, int Nf, int Nfa) {

    double xv[3]= {0.0,1.0,2.0} ;
    for (size_t i = 0; i < N; i++) {
        initParticle(particles+i, Nf, xv);
        particles[i].w = &weights[i];
        particles[i].Nfa = Nfa;
        for (size_t el = 0; el < 2*Nfa; el++) {
            particles[i].xf[el] = i;
        }
        for (size_t el = 0; el < 4*Nfa; el++) {
            particles[i].Pf[el] = 0.0;
        }
    }
}


void cleanup_members(Particle* particles, int N) {
    for(int i = 0; i<N; i++) {
        delParticleMembers(particles+i);
    }
}

int main() {
    // Initialize Input

    //We need to implement benchmark over a range of particle lengths here
    const int N = 100;
    double* weights = NULL;
    
    int N_min = 2*N;
    
    const size_t Nf = 20;
    const size_t Nfa = 10;

    Particle *particles = NULL;
    weights = (double*) malloc(N*sizeof(double));
    particles = (Particle*) malloc(N*sizeof(Particle));
    setup_particles(particles, weights, N, Nf, Nfa);
    
    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&resample_particles)> bench("resample_particles Benchmark");

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    data_loader(particles, N, weights, N_min, 1);
    bench.data_loader = &data_loader;
    bench.add_function(&resample_particles_base, "resample_particles_base", 0.0); // cycles scale exponentially with #Particles!!!
    bench.funcFlops[0] = resample_particles_base_flops(particles, N, weights, N_min, 1);
    bench.funcBytes[0] = 8*resample_particles_base_memory(particles, N, weights, N_min, 1);
    bench.add_function(&resample_particles_dag, "resample_particles_dag", 0.0);
    bench.funcFlops[1] = resample_particles_dag_flops(particles, N, weights, N_min, 1);
    bench.funcBytes[1] = 8*resample_particles_dag_memory(particles, N, weights, N_min, 1);
    //bench.add_function(&resample_particles_orig, "resample_particles_orig", work);
    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(particles, N, weights,N_min, 1);
    // Free memory
    // destroy(A);
    // destroy(x);
    // destroy(y);
    cleanup_members(particles, N);

    Benchmark<decltype(&resample_particles)> bench_scale("resample_particles Benchmark");
    bench_scale.csv_path = "resample_particles_scale_particles.csv";
    bench_scale.csv_output = false;
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    data_loader(particles, N, weights, N_min, 1);
    bench_scale.data_loader = &data_loader;
    bench_scale.add_function(&resample_particles_base, "resample_particles_base", 1.0); 
    bench_scale.add_function(&resample_particles_dag, "resample_particles_dag", 1.0); 
    
    int Np = 100;
    for (int i = 0; i<8; i++) {
        const int Npi = std::pow(2,i)*Np;
        N_min = 2*Npi; // guarantee to do full resample
        double* ws = (double*) aligned_alloc(32, Npi * sizeof(double));
        Particle* ps = (Particle*) aligned_alloc(32, Npi * sizeof(Particle));
        setup_particles(ps, ws, Npi, Nf, Nfa);
        bench_scale.funcFlops[0] = resample_particles_base_flops(ps, Npi, ws, N_min, 1);
        bench_scale.funcBytes[0] = 8*resample_particles_base_memory(ps, Npi, ws, N_min, 1);
        bench_scale.funcFlops[1] = resample_particles_dag_flops(ps, Npi, ws, N_min, 1);
        bench_scale.funcBytes[1] = 8*resample_particles_dag_memory(ps, Npi, ws, N_min, 1);
        bench_scale.run_name = std::to_string(Npi);
        bench_scale.run_benchmark(ps, Npi, ws, N_min, 1);
        cleanup_members(ps, Npi);
        free(ps);
        free(ws);
    }

    bench_scale.details();

    bench_scale.write_csv_details();

    free(weights);
    free(particles);

    return 0;
}
