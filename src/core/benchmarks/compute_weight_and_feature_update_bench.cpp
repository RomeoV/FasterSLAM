#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <iostream>
#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>
#include <immintrin.h>

#include "compute_weight_and_feature_update.h"
#include "linalg.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

Particle* particle_preloader(Vector3d xv, Vector2d xf[], Matrix2d Pf[], const size_t Nfa,double* weight) {
    Particle* particle = newParticle(50, xv);
    for (int i = 0; i< Nfa; i++) {
        set_xfi(particle,xf[i] ,i);
        set_Pfi(particle, Pf[i], i);
    }
    particle->w = weight;
    return particle;
}

void data_loader(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R) {
}

void enforce_symmetry_Pf(Matrix2d Pf[], const size_t Nfa) {
    for ( int i = 0; i<Nfa; i++) {
        Pf[i][2] = Pf[i][1];
    }
}



int main() {
    // // Initialize Input
    // // prepare particle

    // srand(0);
    // Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
    // const size_t Nfa = 20;
    // Vector2d xf[Nfa];
    // Matrix2d Pf[Nfa];

    // double weight_init[2] = {0.1 , 0.0};

    // Vector2d* z = static_cast<Vector2d *>(aligned_alloc(32, Nfa * sizeof(Vector2d)));

    

    // fill_rand(*xf, 2*Nfa, -2.0,2.0);
    // fill_rand(*Pf, 4*Nfa, 0.0,0.2); //Ignore symmetry for computation for now. Pos. semi-definite

    // fill_rand(*z, 2*Nfa, -2.0,2.0);
    // enforce_symmetry_Pf(Pf, Nfa);

    // Particle* particle = particle_preloader(xv,xf,Pf, Nfa, weight_init);

    // const size_t Nfz = 10; // Nfz <= Nfa
    // int idf[Nfz] = {0,1,2,3,4,5,6,7,8,9}; //Unique

    // Matrix2d R = {1,0,0,1};
    
    
    
    // // Compare two methods
    // "functional equality"_test = [&] {
    //     auto is_close = [&](double lhs, double rhs) -> bool {return std::abs(lhs - rhs) < 1e-4;};
        
    //     Vector2d zp_base[Nfa];
    //     Matrix23d Hv_base[Nfa];
    //     Matrix2d Hf_base[Nfa];
    //     Matrix2d Sf_base[Nfa];

    //     Particle* base_particle = particle_preloader(xv,xf,Pf, Nfa, weight_init);

    //     compute_weight_and_feature_update_base(base_particle, z, idf, Nfz, R);

    //     //decltype(&compute_weight_and_feature_update) func[1] = {&compute_weight_and_feature_update_active};

        
    //     compute_weight_and_feature_update_active(particle, z, idf, Nfz, R);
    //     expect(is_close(*(base_particle->w), *(particle->w)));
    //     for (int i = 0; i < Nfa; i++) {
    //         for (int j = 0; j<4;j++) {
    //             expect(is_close(base_particle->Pf[4*i + j], particle->Pf[4*i + j]));
    //         }

    //         for (int j = 0; j<2;j++) {
    //             expect(is_close(base_particle->xf[2*i + j], particle->xf[2*i + j]));
    //         }
    //     }
    // };




    // // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    // Benchmark<decltype(&compute_weight_and_feature_update)> bench("cf_update Benchmark");

    // double work = Nfz*10; // best-case in flops
    // bench.data_loader = data_loader;

    // bench.add_function(&compute_weight_and_feature_update_base, "cf_update_base", work);
    // bench.add_function(&compute_weight_and_feature_update_active, "cf_update_active", work);
    // //bench.add_function(&compute_weight_and_feature_update, "compute_weight_and_feature_update", work);
    
    // int idf_4[1] = {5};
    // bench.run_benchmark(particle, z, idf_4, 1, R);

    // bench.run_benchmark(particle, z, idf, Nfz, R);

    // int idf_2[Nfz] = {0,2,4,6,8,10,12,14,16,18};
    // bench.run_benchmark(particle, z, idf_2, Nfz, R);


    // int idf_3[20] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19};
    // bench.run_benchmark(particle, z, idf_3, 20, R);

    // delParticleMembersAndFreePtr(particle);
    // //bench.destructor_output = false;

    // bench.details();

    // free(z);


    return 0;
}
