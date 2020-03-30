#include "stratified_random.h"  // import file to test
#include <math.h>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    "predict_true"_test = [] {
        given("I have the arguments N, di") = [] {

            int N = 5;
            double di[5];

            when("I call stratified_random(N, di) with a fixed seed") = [&] {
                
                stratified_random(N, di);

                then("I get the filled array di as expected") = [=] {
                    double exact_di[5] = {0.100000000000, 0.368037543431, 0.478876585364,
                                            0.756619844752, 0.959688006695};
                    for (int i = 0; i < N; i++) {
                        double error = fabs(di[i]-exact_di[i]);
                        expect(that % error < 1e-12) << i;
                    }
                };
            };
        };
    };
};
