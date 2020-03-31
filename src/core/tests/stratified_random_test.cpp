#include "stratified_random.h"  // import file to test
#include <math.h>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    srand(1994);
    "predict_true"_test = [] {
        given("I have the arguments N, di") = [] {

            const size_t N = 5;
            double di[5];

            when("I call stratified_random(N, di) with a fixed seed") = [&] {
                
                stratified_random(N, di);

                then("I get the filled array di as expected") = [=] {
                    double exact_di[5] = {0.117519159763, 0.230069380919,
                             0.491273727683, 0.753994618894, 0.985284319420};
                    for (int i = 0; i < N; i++) {
                        double error = fabs(di[i]-exact_di[i]);
                        expect(that % error < 1e-12) << i;
                    }
                };
            };
        };
    };
};
