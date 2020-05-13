#include "compute_weight.h"  // import file to test
#include "compute_jacobians.h"
#include <vector>  // used for test input
#include "ut.hpp"  // import functionality and namespaces from single header file
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

int main() {
  "compute_weight_simple"_test = [] {
    given("I have a particle, map features/indices and observation noises") = [] {
      // prepare particle
      Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
      Particle* particle = newParticle(5, xv);
      Vector2d xf[3] = {{1,0.1},{1,0.2},{1,0.3}}; 
      for(int i=0; i<3; i++){ //! 2d means of EKF in cartesian world coordinates
        set_xfi(particle, xf[i], i);
      }
      Matrix2d Pf[3] = {{1,0,1,0},{1,0,1,0},{1,0,1,0}}; //! covariance matrices for EKF in polar robot coordinates
      for(int i=0; i<3; i++){ //! covariance matrices for EKF (=Extended Kalman Filter) in polar robot coordinates
        set_Pfi(particle, Pf[i], i);
      }

      Vector2d z[3] = {{0,0}, {0,0}, {0,0}};  // vector of map features
      size_t N_z = 3; // number of features.
      int idf[3] = {0,0,0};  // vector of map indices
      Matrix2d R = {1,0,0,1};   // matrix of observation noises

      Vector2d zp[3];
      Matrix23d Hv[3];
      Matrix2d Hf[3];
      Matrix2d Sf[3];
      compute_jacobians_base(particle, idf, 3, R, zp, Hv, Hf, Sf);

      when("I add them") = [&] {
        double weight = compute_weight(particle, z, N_z, idf, R, zp, Hv, Hf, Sf);
        then("I get the right particle weight") = [=] {
          auto is_close = [](auto lhs, auto rhs) -> bool {return lhs - rhs < 1e-14 or rhs - lhs < 1e-14;};
          expect(is_close(weight, 0.000745473));
        };
      };
      
      //! Delete Particle
      delParticleMembersAndFreePtr(particle);
    };
  };

};