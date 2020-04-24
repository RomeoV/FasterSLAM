#include "stratified_resample.h"  // import file to test
#include <cmath>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    srand(1994);
    "stratified_resample"_test = [] {
        given("Working stratified_random function and calling with w, N_w, Neff, keep") = [] {

            const size_t N_w = 5;
            double w[N_w] = {2,4,6,8,10}; //unnormalized weights
            double Neff = 10.0; //Just written to this variable
            size_t keep[N_w] = {100,100,100,100,100}; //Just written to this variable

            when("I call stratified_resample(w, N_w, &Neff, keep) with a fixed seed") = [&] {
                
                stratified_resample(w, N_w, &Neff, keep);

                then("I get the the outputs w, Neff and keep I want") = [=] {
                    double exact_Neff = 4.090909090909;
                    double target_w[N_w] =  {2,4,6,8,10};
                    size_t exact_keep[N_w] = {1,2,3,4,4};
                    for (int i = 0; i < N_w; i++) {
                        double error_w = fabs(w[i]-target_w[i]);
                        expect(that % error_w < 1e-12) << i;
                        expect(that % keep[i] == exact_keep[i]) << i;
                    }
                    expect(that % fabs(Neff - exact_Neff) < 1.0e-12) <<"Neff error";
                };
            };
        };
    };
};
