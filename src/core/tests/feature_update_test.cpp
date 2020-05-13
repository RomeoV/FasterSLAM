#include "feature_update.h"  // import file to test

#include <cmath>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
  "update feature"_test = [] {
    given(
        "A particle with three known features and vehicle position "
      ) = [] {
      // prepare particle
      Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
      Particle* particle = newParticle(3, xv);
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

      Vector2d zp[3] __attribute__ ((aligned(32)));
      Matrix23d Hv[3] __attribute__ ((aligned(32)));
      Matrix2d Hf[3] __attribute__ ((aligned(32)));
      Matrix2d Sf[3] __attribute__ ((aligned(32)));
      compute_jacobians_base(particle, idf, 3, R, zp, Hv, Hf, Sf);

      when("I update features of the particle") = [&](){
        feature_update(particle, z, idf, N_z, R, zp, Hv, Hf, Sf);
        auto is_close = [](double lhs, double rhs) { return fabs(lhs-rhs) < 1e-5; };

        double target_xf[3] = {0.533446,-0.366554,1};
        then("I expect that the xf are equal to the arr elements.") = [=] {
            for (int i = 0; i<3;i++) {
              double xfi = particle->xf[i];
              expect(that % is_close(xfi, target_xf[i]) == true);
            }
        };
        
        double target_Pf[3] = {0.486263,-0.513737,0.486263};
        then("I expect that the Pf are equal to the arr elements.") = [=] {
            for (int i = 0; i<3;i++) {
              double Pfi = particle->Pf[i];
              expect(that % is_close(Pfi, target_Pf[i]) == true);
            }
        };
       
      };
      delParticleMembersAndFreePtr(particle);
    };
  };
}
