#include "pi_to_pi.h"  // import file to test
#include <math.h>
#include <stdio.h>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    "pi_to_pi"_test = [] {
        given("I have the argument angle") = [] {

            double angle = 5*M_PI;

            when("I call pi_to_pi(angle)") = [&] {
                
                double computed_angle = pi_to_pi(angle);

                then("I get the clipped angle in the range (-PI, PI]") = [=] {

                    double exact_angle = M_PI;
                    double error_same = fabs(computed_angle - exact_angle);
                    
                    expect(that % -M_PI < computed_angle and computed_angle <= M_PI);
                    expect(that % error_same < 1e-14);
                };
            };
        };
    };

    "pi_to_pi base"_test = [] {
        given("I have the argument angle") = [] {

            double angle = 5*M_PI;

            when("I call pi_to_pi_base(angle)") = [&] {
                
                double computed_angle = pi_to_pi_base(angle);

                then("I get the clipped angle in the range (-PI, PI]") = [=] {

                    double exact_angle = M_PI;
                    double error_same = fabs(computed_angle - exact_angle);
                    
                    expect(that % -M_PI < computed_angle and computed_angle <= M_PI);
                    expect(that % error_same < 1e-14);
                };
            };
        };
    };

    "pi_to_pi_arr"_test = [] {
        given("I have the argument angles, N") = [] {
            const size_t N=5;
            double angles[N] = {-5*M_PI, -2*M_PI, -M_PI, 1.5*M_PI, 5*M_PI};

            when("I call pi_to_pi(angle)") = [&] {
                
                pi_to_pi_arr(angles,N);

                then("I get the clipped angles in the range (-PI, PI]") = [=] {

                    double exact_angles[N] = {M_PI,0.0,+M_PI,-M_PI/2,M_PI};
                    for (int i = 0; i<N; i++){
                        double error_same = fabs(angles[i] - exact_angles[i]);
                        
                        expect(that % -M_PI < angles[i] and angles[i] <= M_PI) << i;
                        expect(that % error_same < 1e-14) << i;
                    }
                };
            };
        };
    };

    "pi_to_pi_fmod"_test = [] {
        given("I have the argument angles, N") = [] {
            const size_t N=5;
            double angles[N] = {-5*M_PI, -2*M_PI, -M_PI, 1.5*M_PI, 5*M_PI};

            when("I call pi_to_pi_fmod(angle)") = [&] {
                for (size_t i = 0; i < N; i++) {
                    angles[i] = pi_to_pi_fmod(angles[i]);
                }

                then("I get the clipped angles in the range (-PI, PI]") = [=] {

                    double exact_angles[N] = {M_PI,0.0,+M_PI,-M_PI/2,M_PI};
                    for (int i = 0; i<N; i++){
                        double error_same = fabs(angles[i] - exact_angles[i]);
                        
                        expect(that % -M_PI < angles[i] and angles[i] <= M_PI) << i;
                        expect(that % error_same < 1e-14) << i;
                    }
                };
            };
        };
    };

    "pi_to_pi while"_test = [] {
        given("I have the argument angle") = [] {

            double angle = 5*M_PI;

            when("I call pi_to_pi_while(angle)") = [&] {
                
                double computed_angle = pi_to_pi_while(angle);

                then("I get the clipped angle in the range (-PI, PI]") = [=] {

                    double exact_angle = M_PI;
                    double error_same = fabs(computed_angle - exact_angle);
                    
                    expect(that % -M_PI < computed_angle and computed_angle <= M_PI);
                    expect(that % error_same < 1e-14);
                };
            };
        };
    };

};
