#include "particle.h"  // import file to test
#include <cmath>
#include <stdio.h>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

    "particle"_test = [] {
        given("I have the total number of features Nf") = [] {
            const size_t Nf = 10;

            when("I initialize the particle with newParticle(Nf)") = [=] {
                particle* p = newParticle(Nf);

                then("I expect his parameters Nfa, Nf to be 0, Nf.") = [=] {
                    expect(that % p->Nfa == 0);
                    expect(that % p->Nf == Nf);
                };

                then("I expect its size to be a fixed value.") = [=] {
                    size_t p_size = 2*4+8+12*8+2*8+5*8; //2 ints + 1 double + the (3+9) = 12 double arr elements + 2 pointers + 5*8 voids
                    expect(that % sizeof(*p) == p_size);
                };

                when("I want to set a weight to p.w then") = [=] {
                    double weight = 0.5;
                    p->w = weight;
                    then("I expect that the weight is set.") = [=] {
                        expect(that % p->w == weight);
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
                delParticle(p);
            };
        };
    };

};