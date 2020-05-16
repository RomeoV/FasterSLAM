#include "feature_update.h"  // import file to test

#include <cmath>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
/*  "update feature"_test = [] {
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

      Vector2d zp[3];
      Matrix23d Hv[3];
      Matrix2d Hf[3];
      Matrix2d Sf[3];
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
*/
  "feature_update"_test = [] {
    given("A particle with three known features and vehicle position ") = [] {
      // prepare particle

      Vector3d xv = {1.293967823315060, -0.054066219251330, -0.012642858479510};
      
      Particle* particle = newParticle(3, xv);
      
      Vector2d xf[2] = { {3.227460886446243, -25.613382543676146},
                         {25.570128848983597, 3.630650683089399} }; // Transposed from MATLAB
      
      Matrix2d Pf[2] = { {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843},
                         {0.013819565896226, -0.026186052088964, -0.026186052088964, 0.189525459865311} };

      for(int i = 0; i < 2; i++){ //! 2d means of EKF in cartesian world coordinates
        set_xfi(particle, xf[i], i);
        set_Pfi(particle, Pf[i], i);
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
      //compute_jacobians_base(particle, idf, 2, R, zp, Hv, Hf, Sf);

      when("I update features of the particle") = [&](){

        feature_update(particle, z, idf, N_z, R, zp, Hv, Hf, Sf);

        auto is_close = [](double lhs, double rhs) { return fabs(lhs-rhs) < 1e-8; };

        Vector2d target_xf[2] = { {3.470171202213126, -25.656742169761873},
                                  {25.523980084381321, 4.170480246835258} };

        then("I expect that the xf are equal to the arr elements.") = [=] {
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    expect(that % is_close(particle->xf[i*2+j], target_xf[i][j]) == true);
                }
            }
        };
        
        Matrix2d target_Pf[3] = {{0.099441292170602, 0.008341518741166, 0.008341518741166, 0.005737523050281},
                               {0.006932154241992, -0.012983119323818, -0.012983119323818, 0.092241962970569}};
        then("I expect that the Pf are equal to the arr elements.") = [=] {
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    expect(that % is_close(particle->Pf[i*4+j], target_Pf[i][j]) == true);
                }
            }
        };
       
      };

      delParticleMembersAndFreePtr(particle);
    };
  };
}
