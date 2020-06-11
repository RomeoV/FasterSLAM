#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <iostream>
#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>
#include <immintrin.h>

#include "compute_jacobians.h"
#include "linalg.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

Particle* particle_preloader(Vector3d xv, Vector2d xf[], Matrix2d Pf[], const size_t Nfa) {
    Particle* particle = newParticle(50, xv);
    for (int i = 0; i< Nfa; i++) {
        set_xfi(particle, xf[i], i);
        set_Pfi(particle, Pf[i], i);
    }
    return particle;
}

void enforce_symmetry_Pf(Matrix2d Pf[], const size_t Nfa) {
    for ( int i = 0; i<Nfa; i++) {
        Pf[i][2] = Pf[i][1];
    }
}

// I will try to add this as smooth as possible to the benchmark, but for now do this
void set_work(Benchmark<decltype(&compute_jacobians)>& bench, Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) {
    bench.funcFlops[0] = compute_jacobians_base_flops(particle, idf, N_z, R, zp, Hv, Hf, Sf);
    bench.funcBytes[0] = 8 * compute_jacobians_base_memory(particle, idf, N_z, R, zp, Hv, Hf, Sf);
    for (int i = 1; i < bench.numFuncs; i++) {
        bench.funcFlops[i] = compute_jacobians_active_flops(particle, idf, N_z, R, zp, Hv, Hf, Sf);
        bench.funcBytes[i] = 8 * compute_jacobians_active_memory(particle, idf, N_z, R, zp, Hv, Hf, Sf);
    }
}


int main() {
    // Initialize Input
    // create 4 particles... just pretend they are all the same!
    const size_t Nfz = 12; // Nfz <= Nfa

    srand(0);
    Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
    const size_t Nfa = 20;
    Vector2d xf[Nfa];
    Matrix2d Pf[Nfa];
    
    fill_rand(*xf, 2*Nfa, -2.0,2.0);
    fill_rand(*Pf, 4*Nfa, 0.0,0.2); //Ignore symmetry for computation for now. Pos. semi-definite
    enforce_symmetry_Pf(Pf, Nfa);

    Particle* particle = particle_preloader(xv,xf,Pf, Nfa);

    // outputs
    Vector2d zp[Nfa] __attribute__((aligned(32)));
    Matrix23d* Hv = static_cast<Matrix23d *>(aligned_alloc(32, Nfa * sizeof(Matrix23d)));
    Matrix2d* Hf = static_cast<Matrix2d *>(aligned_alloc(32, Nfa * sizeof(Matrix2d)));
    Matrix2d* Sf = static_cast<Matrix2d *>(aligned_alloc(32, Nfa * sizeof(Matrix2d)));

    //Input -> Change this as you like
    int idf[Nfz] __attribute__((aligned(32))) = {0,1,2,3,4,5,6,7,8,9,10,11}; //Unique

    __m256d px = _mm256_set_pd(particle->xv[0],particle->xv[0],particle->xv[0],particle->xv[0]);
    __m256d py = _mm256_set_pd(particle->xv[1],particle->xv[1],particle->xv[1],particle->xv[1]);
    __m256d ptheta = _mm256_set_pd(particle->xv[2],particle->xv[2],particle->xv[2],particle->xv[2]);       
    __m256d xfp0p2 = _mm256_set_pd(particle->xf[2*idf[1]],particle->xf[2*idf[1]],particle->xf[2*idf[1]],particle->xf[2*idf[1]]);  
    __m256d xfp1p3 = _mm256_set_pd(particle->xf[2*idf[1]+1],particle->xf[2*idf[1]+1],particle->xf[2*idf[1]+1],particle->xf[2*idf[1]+1]);
    // particle.Pf is a array of matrices
    __m256d Pf0 = _mm256_load_pd(particle->Pf + 4 * idf[0]); // get matrix for first feature for particle #0
    __m256d Pf1 = _mm256_load_pd(particle->Pf + 4 * idf[0]); // get matrix for first feature for particle #1 = #0
    __m256d Pf2 = _mm256_load_pd(particle->Pf + 4 * idf[0]); // get matrix for first feature for particle #2 = #0
    __m256d Pf3 = _mm256_load_pd(particle->Pf + 4 * idf[0]); // get matrix for first feature for particle #3 = #0            
                            
    __m256d R_vec = _mm256_set_pd(1,0,0,1); // we assume R is the same for all 4 particles
    
    // outputs
    __m256d* zp_dist = static_cast<__m256d *>(aligned_alloc(32, Nfa * sizeof(__m256d)));
    __m256d* zp_angle = static_cast<__m256d *>(aligned_alloc(32, Nfa * sizeof(__m256d)));
    __m256d* Hf0v = static_cast<__m256d *>(aligned_alloc(32, Nfa * sizeof(__m256d)));  
    __m256d* Hf1v = static_cast<__m256d *>(aligned_alloc(32, Nfa * sizeof(__m256d)));  
    __m256d* Hf2v = static_cast<__m256d *>(aligned_alloc(32, Nfa * sizeof(__m256d))); 
    __m256d* Hf3v = static_cast<__m256d *>(aligned_alloc(32, Nfa * sizeof(__m256d)));
    __m256d* s_zeros = static_cast<__m256d *>(aligned_alloc(32, Nfa * sizeof(__m256d))); 
    __m256d* s_ones = static_cast<__m256d *>(aligned_alloc(32, Nfa * sizeof(__m256d)));
    __m256d* s_twos = static_cast<__m256d *>(aligned_alloc(32, Nfa * sizeof(__m256d)));
    __m256d* s_threes = static_cast<__m256d *>(aligned_alloc(32, Nfa * sizeof(__m256d)));

    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&compute_jacobians4x_avx)> bench("compute_jacobians Benchmark AVX");

#ifdef __AVX2__
    bench.add_function(&compute_jacobians4x_avx, "compute_jacobians_avx", 0.0);
#endif

    bench.run_benchmark(px, py, ptheta, xfp0p2, xfp1p3, 
                            Pf0, Pf1, Pf2, Pf3, R_vec,
                            zp_dist, zp_angle,
                            Hf0v, Hf1v, Hf2v, Hf3v,
                            s_zeros, s_ones, s_twos, s_threes);

    bench.details();
 
    delParticleMembersAndFreePtr(particle);
    free(Hv);
    free(Hf);
    free(Sf);

    return 0;
}
