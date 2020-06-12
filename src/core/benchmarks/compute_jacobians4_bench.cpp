#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "particle.h"
#include "linalg.h"
#include "trigonometry.h"
#include "fastrand.h"
#include "compute_jacobians.h"
#include "get_observations.h"

#include <immintrin.h>

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

# define Nfa_start 0
extern __inline __m256d __attribute__((__gnu_inline__, __always_inline__, __artificial__))
_mm256_loadme_m128d (double const *__PH, double const *__PL)
{
  return _mm256_insertf128_pd (_mm256_castpd128_pd256 (_mm_loadu_pd (__PL)),
			       _mm_load_pd (__PH), 1);
}

template<typename func>
std::function<void (int *, int)> random_lambda(func rand_func) {
    return [=] (int* recipient, int N)
    {
        for (int i = 0; i<N; i++) {
            recipient[i] = rand_func();
            
        }
    };
}

__m256i xv_load_mask = _mm256_set_epi64x(9,6,3,0);
__m256d R_vec = _mm256_set_pd(1,0,0,1);
template<typename func>
std::function<void (int *, int)> random_vec_lambda(func rand_func) {
    auto mask = _mm256_set1_epi32 (1);
    return [=] (int* recipient, int N)
    {
        for (int i = 0; i<N; i+=4) {
            _mm256_maskstore_epi32(recipient+i,mask, rand_func());
        }
    };
}


inline void func(Particle* particles, int N, int* idf, int Nfa, double *R, Vector2d* zp, Matrix23d* Hv, Matrix2d* Hf, Matrix2d* Sf) {

    for (int i = 0; i<N; i++) {
        compute_jacobians_base(particles+i, idf, Nfa, R, zp+i*Nfa, Hv+i*Nfa, Hf+i*Nfa, Sf+i*Nfa);
    }
}


void data_loader(Particle* particles, int N, int* idf, int Nfa, double *R, Vector2d* zp, Matrix23d* Hv, Matrix2d* Hf, Matrix2d* Sf) {
    srand(0);
}


inline void vec_func(Particle* particles, int N, int* idf, int Nfa, double *R, Vector2d* zp, Matrix23d* Hv, Matrix2d* Hf, Matrix2d* Sf) {

    // for (int i = 0; i<N; i++) {
    //     compute_jacobians(particles+i, idf, Nfa, R, zp, Hv, Hf, Sf);
    // }
    
    
    __m256d zp_dist, zp_angle, Hf0v,  Hf1v,  Hf2v, Hf3v,s_zeros, s_ones, s_twos, s_threes;
    for (size_t i = 0; i < N; i+=4) {
        __m256d px = _mm256_loadu_pd(particles[i].xv);
        __m256d py = _mm256_loadu_pd(particles[i].xv + 4);
        __m256d ptheta = _mm256_loadu_pd(particles[i].xv+8);
        

        //weights_v = _mm256_load_pd(weights+i);

        for (size_t j = 0; j < Nfa; j++) {
            __m256d xfp0p2 = _mm256_loadme_m128d(particles[i+2].xf + 2 * idf[j], particles[i].xf + 2 * idf[j]);
            __m256d xfp1p3 = _mm256_loadme_m128d(particles[i+3].xf + 2 * idf[j], particles[i+1].xf + 2 * idf[j]);

            __m256d Pf0 = _mm256_load_pd(particles[i].Pf + 4 * idf[j]);
            __m256d Pf1 = _mm256_load_pd(particles[i+1].Pf + 4 * idf[j]);
            __m256d Pf2 = _mm256_load_pd(particles[i+2].Pf + 4 * idf[j]);
            __m256d Pf3 = _mm256_load_pd(particles[i+3].Pf + 4 * idf[j]);

            compute_jacobians4x_avx_inline(px, py, ptheta, xfp0p2, xfp1p3, 
                        Pf0, Pf1, Pf2, Pf3, R_vec,
                        &zp_dist, &zp_angle,
                        &Hf0v,  &Hf1v,  &Hf2v, &Hf3v,
                        &s_zeros, &s_ones, &s_twos, &s_threes);
        }
    }
}


void init_particles_contigous(Particle* particle, double* _xv, double* _Pv, double* weight, const size_t N, const size_t Nf) {
    for (int i = 0; i<N; i++) {
        particle[i].xv = _xv+3*i;
        particle[i].Pv = _Pv+9*i;
        weight[i] = 1.0/N;
        particle[i].w = weight+i;
        initParticle_prealloc(particle+i, Nf, _xv+3*i);
        particle[i].Nfa = Nf;
        fill_rand(particle[i].xf, 2*Nf, -1000.0,1000.0);
        fill_rand(particle[i].Pf, 4*Nf, 0.0,0.9);
    } 
}

void cleanup_members(Particle* particles, int N) {
    for(int i = 0; i<N; i++) {
        delParticleMembers_prealloc(particles+i);
    }
}

int main() {
    const int Nf = 20;
    double lm[2*Nf];
    const int N = 4; // Cannot set to NPARTICLES...

    double xtrue[3] = {0,0,0};
    int N_waypoints = 3;
    double wp[6] = {0,0,1,1,2,2};
    
    
    srand(0);

    srand(0);
    size_t Nf_visible;
    Vector2d *z;
    Vector2d *zf;
    Vector2d *zn;
    int *ftag_visible, *ftag, *da_table;
    double xvs[3*N] __attribute__ ((aligned(32)));
    double Pvs[9*N] __attribute__ ((aligned(32)));
    fill_rand(xvs, 3*N, -0.1,0.1);
    Particle particles[N]  __attribute__ ((aligned(32)));
    double weights[N]  __attribute__ ((aligned(32)));

    const size_t Nfa = 100;

    // Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
    // const size_t Nfa = 20;
    // Vector2d xf[Nfa];
    // Matrix2d Pf[Nfa];

    init_particles_contigous(particles, xvs, Pvs, weights, N, Nfa);

    srand(0);

    Matrix2d R __attribute__((aligned(32))) = {1,0,0,1};

    // outputs
    Vector2d zp[N*Nfa] __attribute__ ((aligned(32)));
    Matrix23d Hv[N*Nfa] __attribute__ ((aligned(32)));
    Matrix2d Hf[N*Nfa] __attribute__ ((aligned(32)));
    Matrix2d Sf[N*Nfa] __attribute__ ((aligned(32)));

    //Input -> Change this as you like
    const size_t Nfz = 20; // Nfz <= Nfa
    int idf[Nfz] __attribute__((aligned(32))) = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19}; //Unique

    Benchmark<decltype(&func)> bench("compute_jacobians Benchmark");

    bench.add_function(&func, "compute_jacobians_base", 0.0);
    bench.add_function(&vec_func, "compute_jacobians_fast", 0.0);

    bench.run_benchmark(particles, N, idf, Nfz, R, zp, Hv, Hf, Sf);
    
}