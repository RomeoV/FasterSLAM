#include "compute_jacobians.h"  // import file to test

#include <vector>  // used for test input
#include "ut.hpp"  // import functionality and namespaces from single header file
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`
#include <cmath>
#include <vector>
#include <functional>

int main() {
    std::vector<std::pair<std::string, std::function<void (Particle*, int[], size_t, Matrix2d, Vector2d[],
                                                           Matrix23d[], Matrix2d[], Matrix2d[])>>>
    compute_jacobians_functions = {
      {"compute_jacobians", compute_jacobians},
      {"compute_jacobians_base", compute_jacobians_base},
      {"compute_jacobians_fast", compute_jacobians_fast},
      //{"compute_jacobians_active", compute_jacobians_active},
      {"compute_jacobians_basic_optimizations", compute_jacobians_basic_optimizations},
#ifdef __AVX2__
      {"compute_jacobians_advanced_optimizations", compute_jacobians_advanced_optimizations},
#endif
      //{"compute_jacobians_simd", compute_jacobians_simd},
      {"compute_jacobians_nik", compute_jacobians_nik},
      {"compute_jacobians_scalar_replacement", compute_jacobians_scalar_replacement},
      {"compute_jacobians_linalg_inplace", compute_jacobians_linalg_inplace}
    };
            
  "compute_jabobians"_test = [](auto NamedFunction) {
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
        //std::cout << "returned " << Hv[0][0] << ", " << Hv[0][1] << ", " << Hv[0][2] << ", " << Hv[0][3] << ", " << Hv[0][4] << ", " << Hv[0][5] << std::endl;

        then("I get the right values back") = [=] {
            auto is_close = [](auto lhs, auto rhs) -> bool {return fabs(lhs - rhs) < 1e-5;};

            Vector2d target_zp[3] = {{1.00499,0.0996687},
                                     {1.00499,0.0996687},
                                     {1.00499,0.0996687}};
            for(int j=0; j<3; j++){
                for(int i=0; i<2; i++){
                    //std::cout << "zp " << zp[j][i] << " target_zp " << target_zp[j][i] << std::endl;
                    expect(is_close(zp[j][i], target_zp[j][i]));
                }
            }
            Matrix23d target_Hv[3] = {{-0.995037,-0.0995037,0,0.0990099,-0.990099,-1},
                                      {-0.995037,-0.0995037,0,0.0990099,-0.990099,-1},
                                      {-0.995037,-0.0995037,0,0.0990099,-0.990099,-1}};
            for(int j=0; j<3; j++){
                for(int i=0; i<6; i++){
                    //std::cout << "Hv " << Hv[j][i] << " target_Hv " << target_Hv[j][i] << std::endl;
                    expect(is_close(Hv[j][i], target_Hv[j][i]));
                }
            }
            Matrix2d target_Hf[4*3] = {{0.995037,0.0995037,-0.0990099,0.990099},
                                       {0.995037,0.0995037,-0.0990099,0.990099},
                                       {0.995037,0.0995037,-0.0990099,0.990099}};
            for(int j=0; j<3; j++){
                for(int i=0; i<4; i++){
                    std::cout << "Hf " << Hf[j][i] << " target_Hf " << target_Hf[j][i] << std::endl;
                    expect(is_close(Hf[j][i], target_Hf[j][i]));
                }
            }
            Matrix2d target_Sf[3] = {{2.08911,-0.10837,0.886667,0.911773},
                                     {2.08911,-0.10837,0.886667,0.911773},
                                     {2.08911,-0.10837,0.886667,0.911773}};
            for(int j=0; j<3; j++){
                for(int i=0; i<4; i++){
                    //std::cout << "Sf " << Sf[j][i] << " target_Sf " << target_Sf[j][i] << std::endl;
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