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
    // prepare particle
    srand(0);
    Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
    const size_t Nfa = 20;
    Vector2d xf[Nfa];
    Matrix2d Pf[Nfa];
    
    fill_rand(*xf, 2*Nfa, -2.0,2.0);
    fill_rand(*Pf, 4*Nfa, 0.0,0.2); //Ignore symmetry for computation for now. Pos. semi-definite
    enforce_symmetry_Pf(Pf, Nfa);

    Particle* particle = particle_preloader(xv,xf,Pf, Nfa);

    Matrix2d R __attribute__((aligned(32))) = {1,0,0,1};

    // outputs
    Vector2d zp[Nfa] __attribute__((aligned(32)));
    Matrix23d* Hv = static_cast<Matrix23d *>(aligned_alloc(32, Nfa * sizeof(Matrix23d)));
    Matrix2d* Hf = static_cast<Matrix2d *>(aligned_alloc(32, Nfa * sizeof(Matrix2d)));
    Matrix2d* Sf = static_cast<Matrix2d *>(aligned_alloc(32, Nfa * sizeof(Matrix2d)));

    //Input -> Change this as you like
    const size_t Nfz = 12; // Nfz <= Nfa
    int idf[Nfz] __attribute__((aligned(32))) = {0,1,2,3,4,5,6,7,8,9,10,11}; //Unique

    // Compare two methods
    "functional equality"_test = [&] {
        auto is_close = [&](double lhs, double rhs) -> bool {return std::abs(lhs - rhs) < 1e-4;};
        
        Vector2d zp_base[Nfa] __attribute__((aligned(32)));
        Matrix23d Hv_base[Nfa] __attribute__((aligned(32)));
        Matrix2d Hf_base[Nfa] __attribute__((aligned(32)));
        Matrix2d Sf_base[Nfa] __attribute__((aligned(32)));

        compute_jacobians_base(particle, idf, Nfz, R, zp_base, Hv_base, Hf_base, Sf_base);
        compute_jacobians_fast(particle, idf, Nfz, R, zp, Hv, Hf, Sf);

        for (int i = 0; i < Nfz; i++) {
            for (int j = 0; j<4;j++) {
                expect(is_close(Hf_base[i][j], Hf[i][j])) <<i << "Hf" <<j;
                expect(is_close(Sf_base[i][j], Sf[i][j])) <<i <<"Sf"<<j << Sf[i][j]<< Sf_base[i][j];
            }

            for (int j = 0; j<6;j++) {
                expect(is_close(Hv_base[i][j], Hv[i][j])) <<i << "Hv"<<j;
            }

            for (int j = 0; j<2;j++) {
                expect(is_close(zp_base[i][j], zp[i][j])) <<i << "zp"<<j << zp[i][j]<< zp_base[i][j];
            }
        }
    };

    auto ymm0 = _mm256_set_pd(3,2,1,0);


    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&compute_jacobians)> bench("compute_jacobians Benchmark");

    bench.add_function(&compute_jacobians_base, "compute_jacobians_base", 0.0);
    bench.add_function(&compute_jacobians_fast, "compute_jacobians_fast", 0.0);
#ifdef __AVX2__
    bench.add_function(&compute_jacobians_advanced_optimizations, "compute_jacobians_jonathan", 0.0);
#endif
    //bench.add_function(&compute_jacobians, "compute_jacobians", work);
    bench.add_function(&compute_jacobians_simd, "compute_jacobians_simd", 0.0);
    bench.add_function(&compute_jacobians_nik, "compute_jacobians_nik", 0.0);
    bench.add_function(&compute_jacobians_scalar_replacement, "compute_jacobians_scalar_replacement", 0.0);
    bench.add_function(&compute_jacobians_linalg_inplace, "compute_jacobians_linalg_inplace", 0.0);

    int idf_4[1] = {5};
    set_work(bench, particle, idf_4, 1, R, zp, Hv, Hf, Sf);
    bench.run_benchmark(particle, idf_4, 1, R, zp, Hv, Hf, Sf);

    set_work(bench, particle, idf, Nfz, R, zp, Hv, Hf, Sf);
    bench.run_benchmark(particle, idf, Nfz, R, zp, Hv, Hf, Sf);

    int idf_2[Nfz] __attribute__((aligned(32))) = {0,2,4,6,8,10,12,14,16,17,18,19};
    set_work(bench, particle, idf_2, Nfz, R, zp, Hv, Hf, Sf);
    bench.run_benchmark(particle, idf_2, Nfz, R, zp, Hv, Hf, Sf);


    int idf_3[20] __attribute__((aligned(32))) = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    set_work(bench, particle, idf_3, 20, R, zp, Hv, Hf, Sf);
    bench.run_benchmark(particle, idf_3, 20, R, zp, Hv, Hf, Sf);


    delParticleMembersAndFreePtr(particle);
    free(Hv);
    free(Hf);
    free(Sf);

    bench.details();

    return 0;
}
