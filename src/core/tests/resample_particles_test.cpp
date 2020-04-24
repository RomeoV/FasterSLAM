#include "resample_particles.h"  // import file to test
#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;


//PLease rewrite this test. resampling (base) is now guaranteed to work.
int main() {
    "resample"_test = [] {
        given("I have three particles with weights [2./3, 1./3, 0.0]") = [] {
            double weights[3] = {2./3, 1./3, 0.};
            const size_t Nf = 5;

            Vector3d zeros = {0,0,0};
            Particle particles[3];
            particles[0].w = &weights[0]; std::copy(particles[0].xv, particles[0].xv+3, zeros); particles[0].Nf = Nf;
            particles[1].w = &weights[1]; std::copy(particles[1].xv, particles[1].xv+3, zeros); particles[1].Nf = Nf;
            particles[2].w = &weights[2]; std::copy(particles[2].xv, particles[2].xv+3, zeros); particles[2].Nf = Nf;
            
            
            for (size_t i = 0; i < 3; i++) {
                double xv_initial[3] =  {0,0,0};
                initParticle(&particles[i], 5,xv_initial);
                particles[i].Nfa = 3;
                for (size_t el = 0; el < 2*3; el++) {
                    particles[i].xf[el] = i;
                }
            }
            when("I resample the particles") = [&] {
                resample_particles(particles, 3, weights, 2, 1);
                
                then("The 0th particle will stay the same, "
                     "the 1st particle will deep copy the 0th "
                     "and the 2nd will deep copy the 1st,"
                     "with the execution order [2, 1, 0]") = [&] {
                    "0th particle"_test = [&]{
                        "xv"_test = [&](const auto& i) {
                            expect(particles[0].xv[i] == 0);
                        } | std::vector{0,1,2};
                        "xf"_test = [&](const auto& i) {
                            expect(particles[0].xf[i] == 0);
                        } | std::vector{0,1,2,3,4,5};
                    };
                    "1st particle"_test = [&]{
                        "xv"_test = [&](const auto& i) {
                            expect(particles[1].xv[i] == 0);
                        } | std::vector{0,1,2};
                        "xf"_test = [&](const auto& i) {
                            expect(particles[1].xf[i] == 0);
                        } | std::vector{0,1,2,3,4,5};
                        "pointer location"_test = [&] {
                            expect(particles[1].xf != particles[0].xf) << "xf";
                            expect(particles[1].Pf != particles[0].Pf) << "Pf";
                        };

                    };
                    "2nd particle"_test = [&]{
                        "xv"_test = [&](const auto& i) {
                            expect(particles[2].xv[i] == 1);
                        } | std::vector{0,1,2};
                        "xf"_test = [&](const auto& i) {
                            expect(particles[2].xf[i] == 1);
                        } | std::vector{0,1,2,3,4,5};
                        "pointer location"_test = [&] {
                            expect(particles[2].xf != particles[1].xf) << "xf";
                            expect(particles[2].Pf != particles[1].Pf) << "Pf";
                        };
                    };
                };
            };
            for (size_t i = 0; i < 3; i++) { delParticleMembers(&particles[i]); }
        };
    };
}
