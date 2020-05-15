#include "compute_jacobians.h"  // import file to test

#include <vector>  // used for test input
#include "ut.hpp"  // import functionality and namespaces from single header file
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`
#include <cmath>

int main() {

  "jabobian_simple"_test = [] {
    given("I have a particle, features and a covariance matrix of observation") = [] {
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
      Matrix2d R = {1,0,0,1};
      when("I compute the jacobians") = [&] {
        Vector2d zp[3] = {{0,0}, {0,0}, {0,0}}; // measurement (range, bearing)
        Matrix23d Hv[3] = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}}; // jacobians of function h (deriv of h wrt pose)
        Matrix2d Hf[3] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}}; // jacobians of function h (deriv of h wrt mean)
        Matrix2d Sf[3] = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}}; // Measurement covariance of feature observation given the vehicle.
        compute_jacobians(particle, idf, N_z, R, // in
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
  };


  "compute_jacobians_test"_test = [] {
      given("I have the arguments x, P, v, R, H") = [] {

          Vector3d xv = {1.293967823315060, -0.054066219251330, -0.012642858479510};
          
          Particle* particle = newParticle(5, xv);
          
          //Vector2d xf[2] = { {3.227460886446243, 25.57012884898359}, 
          //                   {-25.613382543676146, 3.630650683089399} }; 
          Vector2d xf[2] = { {3.227460886446243, -25.613382543676146}, 
                             {25.57012884898359, 3.630650683089399} }; 
         
          Matrix2d Pf[2] = { {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843}, 
                             {0.013819565896226, -0.026186052088964, -0.026186052088964, 0.189525459865311} };

          for(int i = 0; i < 3; i++){
              set_xfi(particle, xf[i], i);
              set_Pfi(particle, Pf[i], i);
          }

          int idf[2] = {0, 1};

          int N_z = 2;
          
          Matrix2d R = {0.010000000000000, 0, 0, 0.000304617419787};

          when("I call compute_jacobians()") = [&] {
              Vector2d zp[2] = {};
              Matrix23d Hv[2] = {};
              Matrix2d Hf[2] = {};
              Matrix2d Sf[2] = {};

              compute_jacobians(/* in -> */ particle, idf, N_z, R, 
                                /* out -> */ zp, Hv, Hf, Sf);

              then("I get the updated values of x and P I want") = [=] {

                  //Vector2d target_zp[2]= { {25.632343755442758, 24.554208046576935},
                  //                         {-1.482649980501718, 0.163276449965079} };
                  Vector2d target_zp[2]= { {25.632343755442758, -1.482649980501718},
                                           {24.554208046576935, 0.163276449965079} };

                  Matrix23d target_Hv[2] = { {-0.075431770172036, 0.997150965525639, 0, 
                                              -0.038902059641499, -0.002942835461779, -1.000000000000000},
                                             {-0.988676196748803, -0.150064579372757, 0,
                                              0.006111562591964, -0.040265041123435, -1.000000000000000} };

                  Matrix2d target_Hf[2] = { {0.075431770172036, -0.997150965525639,
                                             0.038902059641499, 0.002942835461779}, 
                                            {0.988676196748803, 0.150064579372757,
                                             -0.006111562591964, 0.040265041123435} };

                  Matrix2d target_Sf[2] = { {0.020121592762187, -0.000188340825409,
                                             -0.000188340825409, 0.000611567807302},
                                            {0.020006151561333, 0.000043250788361,
                                             0.000043250788361, 0.000625294058576} }; 


                  auto is_close = [](auto lhs, auto rhs) -> bool {return fabs(lhs - rhs) < 1e-8;};

                  for(int i = 0; i < 2; i++){
                      for(int j = 0; j < 2; j++){
                          expect( is_close(zp[i][j], target_zp[i][j]) ) << zp[i][j] << " != " << target_zp[i][j];
                      }
                  }
                  
                  for(int i = 0; i < 2; i++){
                      for(int j = 0; j < 6; j++){
                          expect( is_close(Hv[i][j], target_Hv[i][j]) ) << Hv[i][j] << " != " << target_Hv[i][j];
                      }
                  }
                  
                  for(int i = 0; i < 2; i++){
                      for(int j = 0; j < 4; j++){
                          expect( is_close(Hf[i][j], target_Hf[i][j]) ) << Hf[i][j] << " != " << target_Hf[i][j];
                      }
                  }

                  for(int i = 0; i < 2; i++){
                      for(int j = 0; j < 4; j++){
                          expect(is_close(Sf[i][j], target_Sf[i][j])) << Sf[i][j] << " != " << target_Sf[i][j];
                      }
                  }
                  
                  //! Delete Particle
                  delParticleMembersAndFreePtr(particle);

              };
          };
      };
  };

};
