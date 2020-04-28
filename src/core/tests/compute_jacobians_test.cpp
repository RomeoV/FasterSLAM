#include "compute_jacobians.h"  // import file to test

#include <vector>  // used for test input
#include "ut.hpp"  // import functionality and namespaces from single header file
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

/*! Thrun03 Eq. 38, 39 and 60
 *  http://robots.stanford.edu/papers/Thrun03g.pdf
 * 
 *  Note that the sensor model is defined as follows:
 *  $ g(\theta, s_t) = \begin{bmatrix} \sqrt{(x_{\theta} - x_v)^2 + (y_{\theta} - y_v)^2}  \\
 *                                     \arctan(\frac{y_{\theta} - y_v}{x_{\theta} - x_v}) \end{bmatrix} $
 * 
 *  First computes all predicted observations in relative coordinates to the vehicle.
 *  Then computes the jacobians given a particle state and predict observations. [Compute-Intensive, Switch to mask]
 *  @param[in]   Particle   Particle for which the jacobian should be computed.
 *  @param[in]   idf        Feature indices.
 *  @param[in]   N_z        Number of features.
 *  @param[in]   R          Covariance matrix of observation (diagonal).
 *  @param[out]  zp         vector of predicted observation (given the new vehicle state)
 *  @param[out]  Hv         Jacobian of h wrt vehicle states
 *  @param[out]  Hf         Jacobian of h wrt feature states
 *  @param[out]  Sf         Measurement covariance of feature observation given the vehicle.
 */

int main() {
  "jabobian_simple"_test = [] {
    given("I have a particle, features and a covariance matrix of observation") = [] {
      Vector3d xv = {0,0,0};  //! robot pose: x,y,theta (heading dir)
      Particle* particle = newParticle(5, xv);
      Vector2d xf[5] = {{0,0},{0,0},{0,0},{0,0},{0,0}}; //! 2d means of EKF in cartesian world coordinates
      Matrix2d Pf[5] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}}; //! covariance matrices for EKF in polar robot coordinates
      for(int i=0; i<5; i++){
        set_Pfi(particle, Pf[i], i);
      }
      int idf[5] = {0,0,0,0,0};
      int N_z = 5;
      Matrix2d R = {1,0,0,1};
      when("I compute the jacobians") = [&] {
        Vector2d zp; // measurement (range, bearing)
        Matrix23d Hv; // jacobians of function h (deriv of h wrt pose)
        Matrix2d Hf; // jacobians of function h (deriv of h wrt mean)
        Matrix2d Sf; // Measurement covariance of feature observation given the vehicle.
        compute_jacobians(particle, idf, N_z, R, // in
                &zp, &Hv, &Hf, &Sf // out
                );
        then("I get the right values back") = [=] {
            Vector2d target_zp = {0,0};
            for(int i=0; i<2; i++){
                expect(zp[i] == target_zp[i]);
            }
            Matrix23d target_Hv = {0,0,0,0,0,0};
            for(int i=0; i<6; i++){
                expect(Hv[i] == target_Hv[i]);
            }
            Matrix2d target_Hf = {0,0,0,0};
            for(int i=0; i<4; i++){
                expect(Hf[i] == target_Hf[i]);
            }
            Matrix2d target_Sf = {0,0,0,0};
            for(int i=0; i<4; i++){
                expect(Sf[i] == target_Sf[i]);
            }
        };
      };
    };
  };

  "memory leak"_test = [] {
    // double* foo = (double*)malloc(4*sizeof(double));
    // We don't want this anymore because of the CI pipeline
  };
};
