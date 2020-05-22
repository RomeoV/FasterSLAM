#include "rdtsc_benchmark.h"
#include <iostream>
#include "feature_update.h"
#include "particle.h"
#include "linalg.h"

#include "ut.hpp"
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

void data_loader(Particle *p, Vector2d *z, int *idf, size_t N_idf, double *R, Vector2d *zp, Matrix23d *Hv, Matrix2d *Hf, Matrix2d *Sf) {
    
    Vector3d xv = {1.293967823315060, -0.054066219251330, -0.012642858479510};
    
	copy(xv, 3, p->xv);
    
    Vector2d xf[2] = { {3.227460886446243, -25.613382543676146},
                       {25.570128848983597, 3.630650683089399} }; // Transposed from MATLAB

    Matrix2d Pf[2] = { {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843},
                       {0.013819565896226, -0.026186052088964, -0.026186052088964, 0.189525459865311} };

    for(int i = 0; i < 2; i++){ //! 2d means of EKF in cartesian world coordinates
        set_xfi(p, xf[i], i);
        set_Pfi(p, Pf[i], i);
    }
}

int main() {

    Vector3d xv = {1.293967823315060, -0.054066219251330, -0.012642858479510};

    Particle* p = newParticle(3, xv);

    Vector2d xf[2] = { {3.227460886446243, -25.613382543676146},
                       {25.570128848983597, 3.630650683089399} }; // Transposed from MATLAB

    Matrix2d Pf[2] = { {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843},
                       {0.013819565896226, -0.026186052088964, -0.026186052088964, 0.189525459865311} };

    for(int i = 0; i < 2; i++){ //! 2d means of EKF in cartesian world coordinates
        set_xfi(p, xf[i], i);
        set_Pfi(p, Pf[i], i);
    }

    Vector2d z[2] = { {25.761106705273054, -1.462835729968151},
                      {24.622003552182658, 0.206077227346340} }; // Transposed from MATLAB
        
    size_t N_z = 2; 
    int idf[2] = {0, 1};
    Matrix2d R = {0.010000000000000, 0, 0, 0.000304617419787}; 

    Vector2d zp[2] = { };
    Matrix23d Hv[2] = { };
    Matrix2d Hf[2] = { };
    Matrix2d Sf[2] = { };

    feature_update(p, z, idf, N_z, R, zp, Hv, Hf, Sf);
    
    auto is_close = [](double lhs, double rhs) { return fabs(lhs-rhs) < 1e-8; };

    Vector2d target_xf[2] = { {3.470171202213126, -25.656742169761873},
                              {25.523980084381321, 4.170480246835258} };

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            expect(that % is_close(p->xf[i*2+j], target_xf[i][j]) == true);
        }
    }

    Matrix2d target_Pf[3] = {{0.099441292170602, 0.008341518741166, 0.008341518741166, 0.005737523050281},
                             {0.006932154241992, -0.012983119323818, -0.012983119323818, 0.092241962970569}};

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            expect(that % is_close(p->Pf[i*4+j], target_Pf[i][j]) == true);
        }
    }

    //
    // Benchmark
    //
    // Initialize the benchmark struct by declaring the type of the function you want to benchmark
    Benchmark<decltype(&feature_update)> bench("feature_update Benchmark");
    
    data_loader(&p, z, idf, N_z, R, zp, Hv, Hf, Sf);
    bench.data_loader = data_loader;
    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&feature_update_base, "feature_update_base", 0.0);
    bench.funcFlops[0] = feature_update_base_flops(&p, z, idf, N_z, R, zp, Hv, Hf, Sf);
    bench.funcBytes[0] = feature_update_base_memory(&p, z, idf, N_z, R, zp, Hv, Hf, Sf);
    bench.add_function(&feature_update, "feature_update", 0.0);
    bench.funcFlops[1] = feature_update_active_flops(&p, z, idf, N_z, R, zp, Hv, Hf, Sf);
    bench.funcBytes[1] = feature_update_active_memory(&p, z, idf, N_z, R, zp, Hv, Hf, Sf);

    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    bench.run_benchmark(p, z, idf, N_z, R, zp, Hv, Hf, Sf);


    // Free memory
    delParticleMembersAndFreePtr(p);

    return 0;
}
