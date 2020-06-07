#include "compute_jacobians.h"  // import file to test
#include "linalg.h"
#include "alignment_utils.h"

#include <vector>  // used for test input
#include "ut.hpp"  // import functionality and namespaces from single header file
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`
#include <cmath>
#include <vector>
#include <functional>


auto compute_jacobians_fast_4particles_factory = [](auto compute_jacobians_fast_4particles_func) {
  return [=](
      Particle* particle, int idf[],
      size_t N_z, Matrix2d R,
      Vector2d zp[], Matrix23d Hv[],
      Matrix2d Hf[], Matrix2d Sf[]) {
      Particle* particle4[4] = {particle, particle, particle, particle};
      Vector2d*  zp4[4] = {(Vector2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Vector2d ), 32)),
                           (Vector2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Vector2d ), 32)),
                           (Vector2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Vector2d ), 32)),
                           (Vector2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Vector2d ), 32))}; // measurement (range, bearing)
      Matrix23d* Hv4[4] = {(Matrix23d*) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix23d), 32)), 
                           (Matrix23d*) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix23d), 32)), 
                           (Matrix23d*) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix23d), 32)), 
                           (Matrix23d*) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix23d), 32))}; // jacobians of function h (deriv of h wrt pose)
      Matrix2d*  Hf4[4] = {(Matrix2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix2d ), 32)),
                           (Matrix2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix2d ), 32)),
                           (Matrix2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix2d ), 32)),
                           (Matrix2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix2d ), 32))}; // jacobians of function h (deriv of h wrt mean)
      Matrix2d*  Sf4[4] = {(Matrix2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix2d ), 32)),
                           (Matrix2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix2d ), 32)),
                           (Matrix2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix2d ), 32)),
                           (Matrix2d* ) aligned_alloc(32, aligned_alloc_size(N_z*sizeof(Matrix2d ), 32))}; // Measurement covariance of feature observation given the vehicle.

      compute_jacobians_fast_4particles_func(particle4, idf, N_z, R, zp4, Hv4, Hf4, Sf4);

      "compute_jacobians_fast_4particles internal equality"_test = [&] {
        for (size_t p = 1; p < 4; p++) {  // all particles should be the same, compare 0 to 1-3
          for (size_t j = 0; j < N_z; j++) {  // we compute 3 matrices
            for (size_t el = 0; el < 2; el++) {  // compare each element in matrix
              expect(that % fabs(zp4[p][j][el] - zp4[0][j][el]) < 1e-14);
            }
            //expect(that % true == false) << "It fails";
            for (size_t el = 0; el < 6; el++) {
              // Unused, so we don't always calculate it!
              expect(that % fabs(Hv4[p][j][el] - Hv4[0][j][el]) < 1e-14);
            }
            for (size_t el = 0; el < 4; el++) {
              expect(that % fabs(Hf4[p][j][el] - Hf4[0][j][el]) < 1e-14);
            }
            for (size_t el = 0; el < 4; el++) {
              expect(that % fabs(Sf4[p][j][el] - Sf4[0][j][el]) < 1e-14);
            }
          }
        }
      };
      for (size_t j = 0; j < N_z; j++) {  // copy back to output to check
        copy(zp4[0][j], 2, zp[j]);
        copy(Hv4[0][j], 6, Hv[j]);
        copy(Hf4[0][j], 4, Hf[j]);
        copy(Sf4[0][j], 4, Sf[j]);
      }

      for (size_t p = 0; p < 4; p++) {
        free(zp4[p]);
        free(Hv4[p]);
        free(Hf4[p]);
        free(Sf4[p]);
      }
  };
};

int main() {
    std::vector<std::pair<std::string, std::function<void (Particle*, int[], size_t, Matrix2d, Vector2d[],
                                                           Matrix23d[], Matrix2d[], Matrix2d[])>>>
    compute_jacobians_functions = {
      {"compute_jacobians", compute_jacobians},
      {"compute_jacobians_base", compute_jacobians_base},
      {"compute_jacobians_fast", compute_jacobians_fast},
      {"compute_jacobians_fast_4particles", compute_jacobians_fast_4particles_factory(compute_jacobians_fast_4particles)},
      {"compute_jacobians_fast_4particles_fullavx", compute_jacobians_fast_4particles_factory(compute_jacobians_fast_4particles_fullavx)},
      //{"compute_jacobians_active", compute_jacobians_active},
      {"compute_jacobians_basic_optimizations", compute_jacobians_basic_optimizations},
// #ifdef __AVX2__ // This test fails for me (Philipp L.)
//       {"compute_jacobians_advanced_optimizations", compute_jacobians_advanced_optimizations},
// #endif
      //{"compute_jacobians_simd", compute_jacobians_simd},
      {"compute_jacobians_nik", compute_jacobians_nik},
      {"compute_jacobians_scalar_replacement", compute_jacobians_scalar_replacement},
      {"compute_jacobians_linalg_inplace", compute_jacobians_linalg_inplace}
    };

            
  "compute_jabobians"_test = [&](auto NamedFunction) {
    given("I have a particle, features and a covariance matrix of observation") = [&] {
      // prepare particle
      Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
      Particle* particle = newParticle(5, xv);
      Vector2d xf[3] = {{1,0.1},
                        {1,0.2},
                        {1,0.3}}; 
      for(int i=0; i<3; i++){ //! 2d means of EKF in cartesian world coordinates
        set_xfi(particle, xf[i], i);
      }
      Matrix2d Pf[3] = {{1,0,1,0},
                        {1,0,1,0},
                        {1,0,1,0}}; //! covariance matrices for EKF in polar robot coordinates
      for(int i=0; i<3; i++){ //! covariance matrices for EKF (=Extended Kalman Filter) in polar robot coordinates
        set_Pfi(particle, Pf[i], i);
      }
      int idf[3] = {0,0,0};
      int N_z = 3;
      Matrix2d R __attribute__((aligned(32))) = {1,0,0,1};
      when("I compute the jacobians using " + NamedFunction.first) = [&] {
        Vector2d zp[3] __attribute__((aligned(32))) = {{0,0}, {0,0}, {0,0}}; // measurement (range, bearing)
        Matrix23d Hv[3] __attribute__((aligned(32))) = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}}; // jacobians of function h (deriv of h wrt pose)
        Matrix2d Hf[3] __attribute__((aligned(32))) = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}}; // jacobians of function h (deriv of h wrt mean)
        Matrix2d Sf[3] __attribute__((aligned(32))) = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}}; // Measurement covariance of feature observation given the vehicle.
        NamedFunction.second(particle, idf, N_z, R, // in
                zp, Hv, Hf, Sf // out
                );

        then("I get the right values back") = [=] {
            auto is_close = [](auto lhs, auto rhs) -> bool {return fabs(lhs - rhs) < 1e-5;};

            Vector2d target_zp[3] = {{1.00499,0.0996687},
                                     {1.00499,0.0996687},
                                     {1.00499,0.0996687}};
            for(int j=0; j<3; j++){
                for(int i=0; i<2; i++){
                    expect(is_close(zp[j][i], target_zp[j][i]));
                }
            }
            Matrix23d target_Hv[3] = {{-0.995037,-0.0995037,0,0.0990099,-0.990099,-1},
                                      {-0.995037,-0.0995037,0,0.0990099,-0.990099,-1},
                                      {-0.995037,-0.0995037,0,0.0990099,-0.990099,-1}};
            for(int j=0; j<3; j++){
                for(int i=0; i<6; i++){
                    // Unused, so we don't always calculate it!
                    expect(is_close(Hv[j][i], target_Hv[j][i]));
                }
            }
            Matrix2d target_Hf[4*3] = {{0.995037,0.0995037,-0.0990099,0.990099},
                                       {0.995037,0.0995037,-0.0990099,0.990099},
                                       {0.995037,0.0995037,-0.0990099,0.990099}};
            for(int j=0; j<3; j++){
                for(int i=0; i<4; i++){
                    expect(is_close(Hf[j][i], target_Hf[j][i]));
                }
            }
            Matrix2d target_Sf[3] = {{2.08911,-0.10837,0.886667,0.911773},
                                     {2.08911,-0.10837,0.886667,0.911773},
                                     {2.08911,-0.10837,0.886667,0.911773}};
            for(int j=0; j<3; j++){
                for(int i=0; i<4; i++){
                    expect(is_close(Sf[j][i], target_Sf[j][i]));
                }
            }
            //! Delete Particle
            delParticleMembersAndFreePtr(particle);
        };
      };
    };
  } | compute_jacobians_functions;
}
