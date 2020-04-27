#include "resample_particles.h"  // import file to test
#include "ut.hpp"
#include <vector>
using namespace boost::ut;
using namespace boost::ut::bdd;


//PLease rewrite this test. resampling (base) is now guaranteed to work.
int main() {
    "resample"_test = [] {
        given("I have three particles with weights [2./3, 1./3, 0.0]") = [] {
            double weights[3] = {2./3, 1./3, 0.};
            const size_t Nf = 5;

            Particle particles[3];

            Vector3d zeros = {0.,0.,0.};
            particles[0].w = &weights[0];
            particles[0].Nf = Nf;

            Vector3d ones = {1.,1.,1.};
            particles[1].w = &weights[1];
            particles[1].Nf = Nf;

            Vector3d twos = {2.,2.,2.};
            particles[2].w = &weights[2];
            particles[2].Nf = Nf;
            
            
            for (size_t i = 0; i < 3; i++) {
                auto xv = std::vector{zeros, ones, twos};
                initParticle(&particles[i], 5, xv[i]);
                particles[i].Nfa = 3;
                for (size_t el = 0; el < 2*3; el++) {
                    particles[i].xf[el] = i;
                }
            }

            when("I resample the particles") = [&] {
                resample_particles(particles, 3, weights, 2, 1);  // Nmin = 2 > \frac{1}{\sum{w_i^2}} = 1/((2/3)^2 + (1/3)^2 + 0) = 9/5
                
                then("The 0th particle will stay the same, "
                     "the 1st particle will deep copy the 0th "
                     "and the 2nd will deep copy the 1st, "
                     "with the execution order [2, 1, 0]") = [&] {
                    "0th particle"_test = [&]{
                        "xv"_test = [&](const auto& i) {
                            expect(that % particles[0].xv[i] == 0);
                        } | std::vector{0,1,2};
                        "xf"_test = [&](const auto& i) {
                            expect(that % particles[0].xf[i] == 0);
                        } | std::vector{0,1,2,3,4,5};
                    };
                    "1st particle"_test = [&]{
                        "xv"_test = [&](const auto& i) {
                            expect(that % particles[1].xv[i] == 0);
                        } | std::vector{0,1,2};
                        "xf"_test = [&](const auto& i) {
                            expect(that % particles[1].xf[i] == 0);
                        } | std::vector{0,1,2,3,4,5};
                        "pointer location"_test = [&] {
                            expect(that % particles[1].xf != particles[0].xf) << "xf";
                            expect(that % particles[1].Pf != particles[0].Pf) << "Pf";
                        };

                    };
                    "2nd particle"_test = [&]{
                        "xv"_test = [&](const auto& i) {
                            expect(that % particles[2].xv[i] == 1);
                        } | std::vector{0,1,2};
                        "xf"_test = [&](const auto& i) {
                            expect(that % particles[2].xf[i] == 1);
                        } | std::vector{0,1,2,3,4,5};
                        "pointer location"_test = [&] {
                            expect(particles[2].xf != particles[1].xf) << "xf";
                            expect(particles[2].Pf != particles[1].Pf) << "Pf";
                        };
                    };
                };
            };
            
            when("I call resample again, but Nmin < 3") = [&] {
                resample_particles(particles, 3, weights, 2, 1);  // Nmin = 2 < \frac{1}{\sum{w_i^2}} = 1/((1/3)^2 + (1/3) ^ 2 + (1/3)^2) = 3

                then("He will not resample again") = [&] {
                    "particle 0"_test = [&] (size_t i) {
                        expect(that % particles[0].xv[i] == 0);
                    } | std::vector{0,1,2};

                    "particle 1"_test = [&] (size_t i) {
                        expect(that % particles[1].xv[i] == 0);
                    } | std::vector{0,1,2};

                    "particle 2"_test = [&] (size_t i) {
                        expect(that % particles[2].xv[i] == 1);
                    } | std::vector{0,1,2};
                };
            };
            for (size_t i = 0; i < 3; i++) { delParticleMembers(&particles[i]); }
        };
    };
}
