#include "predict.h"  // import file to test
#include <math.h>
#include <random>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    "predict"_test = [] {
        given("I have a particle at position (x=1,y=1,theta=pi/2)") = [] {
            double initial_pos[3] = {1, 1, M_PI/2};
            Particle p{.xv = {1, 1, M_PI/2}};

            when("I give control inputs for a rectangle") = [=]() mutable {
                double Q[4] = {0, 0, 0, 0};

                const double V = random()*1.0/RAND_MAX;
                const double dt = random()*1.0/RAND_MAX;
                double intermediate_pos[3] = {1-2*V*dt, 1 - V*dt, -M_PI/2};

                predict(p, 2*V, M_PI/2 + 0*M_PI, Q, dt);
                predict(p, 1*V, M_PI/2 + 2*M_PI, Q, dt);  // +2PI to check invariance to values outside of (-PI,PI]
                then("After half way, I'm in a different place") = [&] {
                    "position"_test = [&](size_t i) {
                        expect(that % fabs(intermediate_pos[i] - p.xv[i]) < 1e-14) << p.xv[i] << "!= " << intermediate_pos[i];
                    } | std::vector{0,1,2};
                };

                predict(p, 2*V, M_PI/2 + 4*M_PI, Q, dt);
                predict(p, 1*V, M_PI/2 - 2*M_PI, Q, dt);
                then("In the end I end up in the same place as before") = [&] {
                    "position"_test = [&](size_t i) {
                        expect(that % fabs(initial_pos[i] - p.xv[i]) < 1e-14) << p.xv[i];
                    } | std::vector{0,1,2};
                };
            };

            when("I give control inputs for a pentagon") = [=]() mutable {
                double Q[4] = {0, 0, 0, 0};

                const double V = random()*1.0/RAND_MAX;
                const double dt = random()*1.0/RAND_MAX;
                predict(p, 1*V, M_PI*2./5 + 0*M_PI, Q, dt);
                predict(p, 1*V, M_PI*2./5 + 6*M_PI, Q, dt);
                predict(p, 1*V, M_PI*2./5 + 2*M_PI, Q, dt);
                then("After 3/5 of the way, I'm in a different place") = [&] {
                    "position"_test = [&](size_t i) {
                        expect(that % fabs(initial_pos[i] - p.xv[i]) > 1e-14) << p.xv[i];
                    } | std::vector{0,1,2};
                };

                predict(p, 1*V, M_PI*2./5 - 8*M_PI, Q, dt);
                predict(p, 1*V, M_PI*2./5 - 2*M_PI, Q, dt);
                then("I end up in the same place as before") = [&] {
                    "position"_test = [&](size_t i) {
                        expect(that % fabs(initial_pos[i] - p.xv[i]) < 1e-14) << p.xv[i];
                    } | std::vector{0,1,2};
                };
            };
        };
    };
}
