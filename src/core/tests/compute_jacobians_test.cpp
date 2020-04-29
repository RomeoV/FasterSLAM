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
      int idf[3] = {0,0,0};
      int N_z = 3;
      Matrix2d R = {1,0,0,1};
      when("I compute the jacobians") = [&] {
        Vector2d zp[2*3] = {0,0, 0,0, 0,0}; // measurement (range, bearing)
        Matrix23d Hv[6*3] = {0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0}; // jacobians of function h (deriv of h wrt pose)
        Matrix2d Hf[4*3] = {0,0,0,0, 0,0,0,0, 0,0,0,0}; // jacobians of function h (deriv of h wrt mean)
        Matrix2d Sf[4*3] = {0,0,0,0, 0,0,0,0, 0,0,0,0}; // Measurement covariance of feature observation given the vehicle.
        compute_jacobians(particle, idf, N_z, R, // in
                zp, Hv, Hf, Sf // out
                );
        
        then("I get the right values back") = [=] {
            auto is_close = [](auto lhs, auto rhs) -> bool {return lhs - rhs < 1e-14 or rhs - lhs < 1e-14;};

            Vector2d target_zp[2*3] = {1.00499,0.0996687,1.00499,0.0996687,1.00499,0.0996687};
            for(int j=0; j<3; j++){
                for(int i=0; i<2; i++){
                    expect(is_close(zp[j][i], target_zp[j][i]));
                }
            }
            Matrix23d target_Hv[6*3] = {-0.995037,-0.0995037,0,0.0990099,-0.990099,-1,-0.995037,-0.0995037,0,0.0990099,-0.990099,-1,-0.995037,-0.0995037,0,0.0990099,-0.990099,-1};
            for(int j=0; j<3; j++){
                for(int i=0; i<6; i++){
                    expect(is_close(Hv[j][i], target_Hv[j][i]));
                }
            }
            Matrix2d target_Hf[4*3] = {0.995037,0.0995037,-0.0990099,0.990099,0.995037,0.0995037,-0.0990099,0.990099,0.995037,0.0995037,-0.0990099,0.990099};
            for(int j=0; j<3; j++){
                for(int i=0; i<4; i++){
                    expect(is_close(Hf[j][i], target_Hf[j][i]));
                }
            }
            Matrix2d target_Sf[4*3] = {2.08911,-0.10837,0.886667,0.911773,2.08911,-0.10837,0.886667,0.911773,2.08911,-0.10837,0.886667,0.911773};
            for(int j=0; j<3; j++){
                for(int i=0; i<4; i++){
                    expect(is_close(Sf[j][i], target_Sf[j][i]));
                }
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
