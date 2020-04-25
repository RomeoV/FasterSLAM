#include "compute_steering.h"  // import file to test
#include <math.h>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    "compute_steering"_test = [] {
        given("I have the arguments x, wp, N_wp, minD, rateG=0.1, maxG, dt, iwp, G") = [] {

            double x[3] = {0,0,0};
            const int N_wp = 3;
            double wp[2*N_wp] = {0,0,1,1,2,2};
            double minD = 0.1;
            double rateG = 0.1;
            double maxG = 1.5;
            double dt = 1.0;
            int iwp = 1;
            double G = 0.0;

            when("I call compute_steering(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G) with a fixed seed") = [&] {
                
                compute_steering(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);

                then("I get the outputs iwp and G I want") = [=] {
                    expect(that % iwp==1);
                    expect(that % fabs(G-0.10000)<1.0e-12);
                };
            };
        };

        given("I have the arguments x, wp, N_wp, minD, rateG=1.0, maxG, dt, iwp, G") = [] {

            double x[3] = {0,0,0};
            const int N_wp = 3;
            double wp[2*N_wp] = {0,0,1,1,2,2};
            double minD = 0.1;
            double rateG = 1.0;
            double maxG = 1.5;
            double dt = 1.0;
            int iwp = 1;
            double G = 0.0;

            when("I call compute_steering(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G) with a fixed seed") = [&] {
                
                compute_steering(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);

                then("I get the outputs iwp and G I want") = [=] {
                    expect(that % iwp==1);
                    expect(that % fabs(G-0.785398163397)<1.0e-11);
                };
            };
        };

        given("I have the arguments x, wp, N_wp, minD, rateG=1.0, maxG=0.5, dt, iwp, G") = [] {

            double x[3] = {0,0,0};
            const int N_wp = 3;
            double wp[2*N_wp] = {0,0,1,1,2,2};
            double minD = 0.1;
            double rateG = 1.0;
            double maxG = 0.5;
            double dt = 1.0;
            int iwp = 1;
            double G = 0.0;

            when("I call compute_steering(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G) with a fixed seed") = [&] {
                
                compute_steering(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);

                then("I get the outputs iwp and G I want") = [=] {
                    expect(that % iwp==1);
                    expect(that % fabs(G-0.5)<1.0e-11);
                };
            };
        };

        given("I have the arguments x, wp, N_wp, minD=2.0, rateG=1.0, maxG=0.5, dt, iwp, G") = [] {

            double x[3] = {0,0,0};
            const int N_wp = 3;
            double wp[2*N_wp] = {0,0,1,1,2,2};
            double minD = 2.5;
            double rateG = 1.0;
            double maxG = 0.5;
            double dt = 1.0;
            int iwp = 1;
            double G = 0.0;

            when("I call compute_steering(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G) with a fixed seed") = [&] {
                
                compute_steering(x, wp, N_wp, minD, rateG, maxG, dt, &iwp, &G);

                then("I get the outputs iwp and G I want") = [=] {
                    expect(that % iwp==2);
                    expect(that % fabs(G-0.5)<1.0e-11);
                };
            };
        };
    };
};