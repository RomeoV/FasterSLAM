#include "particle.h"  // import file to test
#include <cmath>
#include <stdio.h>
#include <iostream>
#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

    "particle"_test = [] {
        given("I have the total number of features Nf") = [] {
            const size_t Nf = 10;

            when("I initialize the particle with newParticle(Nf)") = [=] {
                Particle* p = newParticle(Nf);

                then("I expect his parameters Nfa, Nf to be 0, Nf.") = [=] {
                    expect(that % p->Nfa == 0);
                    expect(that % p->Nf == Nf);
                };

                then("I expect its size to be a fixed value.") = [=] {
                    size_t p_size = 4*4+12*8+3*8+6*8; //(3+1empty) ints + the (3+9) = 12 double arr elements + 3 pointers + 6*8 voids
                    expect(that % sizeof(*p) == p_size);
                };

                when("I want to set a scalar weight to p.w then (You should use an array indexer instead, see test below)") = [=] {
                    double weight = 0.5;
                    (p->w) = &weight;
                    then("I expect that the weight is set.") = [=] {
                        expect(that % *(p->w) == weight);
                    };
                };

                when("I want to set a double arr[3] as p.xv then") = [=] {
                    double arr[3] = {1,2,3};
                    p->set_xv(p, arr);
                    then("I expect that the set parameters are equal to the arr elements.") = [=] {
                        for (int i = 0; i<3;i++) {
                            expect(that % p->xv[i] == arr[i]);
                        }
                    };
                    then("...they should also equal in a different way of setting.") = [=] {
                        double arr[3] = {4,5,6};
                        set_xv(p, arr);
                        for (int i = 0; i<3;i++) {
                            expect(that % p->xv[i] == arr[i]);
                        }
                    };
                };

                when("I want to set a double arr[9] as p.Pv then") = [=] {
                    double arr[9] = {1,2,3,4,5,6,7,8,9};
                    p->set_Pv(p, arr);
                    then("I expect that the set parameters are equal to the arr elements.") = [=] {
                        for (int i = 0; i<3;i++) {
                            expect(that % p->Pv[i] == arr[i]);
                        }
                    };
                    then("...they should also equal in a different way of setting.") = [=] {
                        double arr[9] = {10,11,12,13,14,15,16,17,18};
                        set_Pv(p, arr);
                        for (int i = 0; i<3;i++) {
                            expect(that % p->Pv[i] == arr[i]);
                        }
                    };
                };

                when("I want to set a double arr[2] to p.xf at index k then") = [=] {
                    double arr[2] = {1,2};
                    const size_t k = 2;
                    p->set_xfi(p, arr, k);
                    then("I expect that the set parameters are equal to the arr elements.") = [=] {
                        for (int i = 0; i<2;i++) {
                            expect(that % p->xf[2*k+i] == arr[i]);
                        }
                    };
                    then("...they should also equal in a different way of setting.") = [=] {
                        double arr[2] = {3,4};
                        set_xfi(p, arr, k);
                        for (int i = 0; i<2;i++) {
                            expect(that % p->xf[2*k+i] == arr[i]);
                        }
                    };
                    then("I also might want to retrieve a pointer to that location.") = [=] {
                        double arr[2] = {3,4};
                        set_xfi(p, arr, k);
                        double* xfi = (p->xf) + 2*k;
                        for (int i = 0; i<2;i++) {
                            expect(that % xfi[i] == arr[i]);
                        }
                    };
                    then("I might want to know the number of features that are filled (we usually fill in order).") = [=] {
                        expect(that % p->Nfa == 2);
                    };
                };

                when("I want to set a double arr[4] to p.Pf at index k then") = [=] {
                    double arr[4] = {1,2,3,4};
                    const size_t k = 5;
                    p->set_Pfi(p, arr, k);
                    then("I expect that the set parameters are equal to the arr elements.") = [=] {
                        for (int i = 0; i<4;i++) {
                            expect(that % p->Pf[4*k+i] == arr[i]);
                        }
                    };
                    then("...they should also equal in a different way of setting.") = [=] {
                        double arr[4] = {5,6,7,8};
                        set_Pfi(p, arr, k);
                        for (int i = 0; i<4;i++) {
                            expect(that % p->Pf[4*k+i] == arr[i]);
                        }
                    };
                    then("I also might want to retrieve a pointer to that location.") = [=] {
                        double arr[4] = {5,6,7,8};
                        set_Pfi(p, arr, k);
                        double* Pfi = (p->Pf) + 4*k;
                        for (int i = 0; i<4;i++) {
                            expect(that % Pfi[i] == arr[i]);
                        }
                    };
                    then("I might want to know the number of features that are filled (we usually fill in order).") = [=] {
                        expect(that % p->Nfa == 5);
                    };
                };
                //! Delete Particle
                delParticleMembersAndFreePtr(p);
            };
        };
    };

    "particle_array"_test = [] {
        given("I have the total number of features Nf, the number of particles N") = [] {
            const size_t Nf = 10;
            const size_t N = 5;

            when("I want to initialize the set of particles and an array of corresponding weights") = [=] {
                double weights[N];
                Particle particles[N];

                for (int i = 0; i<N;i++) {
                    weights[i] = 1.0/N;
                    initParticle(particles+i, Nf, i); //Not sure if we really need the index in the particle, but I keep it for now
                    particles[i].w = weights+i;
                }
                then("I want to access the particle weights in an aligned way.") = [=] {
                    for (int i = 0; i<N; i++) {
                        expect(that % weights[i] == particles[i].w[0]);
                    }
                };
                for (int i = 0; i<N; i++) {
                    delParticleMembers(&particles[i]);
                }
            };
            
            when("I want to initialize the set of particles and an array of corresponding weights and update an element in the weights array") = [=] {
                double weights[N];
                Particle particles[N];

                for (int i = 0; i<N;i++) {
                    weights[i] = 1.0/N;
                    initParticle(particles+i, Nf, i);
                    particles[i].w = weights+i;
                }
                //Update
                double new_weight = 0.8;
                weights[2] = new_weight;
                then("I want the update to be seen in the particle weight.") = [=] {
                    expect(that % *(particles[2].w) == new_weight);
                    expect(that % weights[2] == new_weight);
                };
                for (int i = 0; i<N; i++) {
                    delParticleMembers(&particles[i]);
                }
            };

            when("I want to initialize the set of particles and an array of corresponding weights and update the weight in a particle") = [=] {
                double weights[N];
                Particle particles[N];

                for (int i = 0; i<N;i++) {
                    weights[i] = 1.0/N;
                    initParticle(particles+i, Nf, i);
                    particles[i].w = weights+i;
                }
                //Update
                double new_weight = 0.5;
                *(particles[3].w) = new_weight;
                then("I want the update to be seen in the particle weight.") = [=] {
                    expect(that % *(particles[3].w) == new_weight);
                    expect(that % weights[3] == new_weight);
                };
                for (int i = 0; i<N; i++) {
                    delParticleMembers(&particles[i]);
                }
            };
            
        };
    };
};
