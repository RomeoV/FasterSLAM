#include "compute_weight.h"  // import file to test

#include <vector>  // used for test input
#include "ut.hpp"  // import functionality and namespaces from single header file
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

/*! Thrun03 Eq. 61
    Compute particle weight for sampling.
    Uses compute_jacobians.
    @param[out] particle Particle whose weight is calculated
    @param[in]  z        vector of map features, calculated by data_associate_known
    @param[in]  idf      vector of map indices, calculated by data_associate_known, used for jacobians 
    @param[in]  R        matrix of observation noises, metres and radians
 */

int main() {
  "compute_weight_simple"_test = [] {
    given("I have a particle, map features/indices and observation noises") = [] {
      Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
      Particle* particle = newParticle(5, xv);
      Vector2d xf[5] = {{0,0},{0,0},{0,0},{0,0},{0,0}}; //! 2d means of EKF in cartesian world coordinates
      Matrix2d Pf[5] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}; //! covariance matrices for EKF in polar robot coordinates
      for(int i=0; i<5; i++){
        set_Pfi(particle, Pf[i], i);
      }

      Vector2d z[3] = {{0,0}, {0,0}, {0,0}};  // vector of map features
      size_t N_z = 3; // number of features.
      int idf[3] = {0,0,0};  // vector of map indices
      Matrix2d R = {1,0,0,1};   // matrix of observation noises

      when("I add them") = [&] {
        compute_weight(particle, z, N_z, idf, R);
        then("I get the right particle weight") = [=] {
          double* weight = particle->w; // todo: why is w a pointer?
          expect(weight[0] == 0);
        };
      };
    };
  };

  "memory leak"_test = [] {
    // double* foo = (double*)malloc(4*sizeof(double));
    // We don't want this anymore because of the CI pipeline
  };
};