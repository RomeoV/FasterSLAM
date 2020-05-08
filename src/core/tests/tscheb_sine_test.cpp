#include "tscheb_sine.h"  // import file to test
#include "pi_to_pi.h"
#include <cmath>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    "tscheb_sine"_test = [] {
        given("I have 100 random angles in [-10PI, 10PI]") = [] {
            double angles[100];
            auto rand_angle = []{ return rand()/RAND_MAX * 4. * M_PI - 2. * M_PI; };
            for (size_t i = 0; i < 100; i++) { angles[i] = rand_angle(); }

            when("I calculate the cmath::sin, the tscheb_fsine and the tscheb_dsine") = [&] {
                double sin_results[100];
                double fsine_results[100];
                double dsine_results[100];
                double dsine_vec_results[100];
                double dsine_unrolled_results[100];

                for (size_t i = 0; i < 100; i++) {
                    sin_results[i] = sin(angles[i]);
                    fsine_results[i] = tscheb_fsine(angles[i], false);
                    dsine_results[i] = tscheb_dsine(angles[i], false);
                }

                // normalize angles for methods that can only deal with normalized angles
                for (size_t i = 0; i < 100; i++) {
                    angles[i] = pi_to_pi(angles[i]);
                }
                tscheb_dsines(angles, 100, dsine_vec_results);
                tscheb_dsines_unrolled(angles, 100, dsine_unrolled_results);

                then("I expect each approximation to be close to cmath::sin, up to some tolerance") = [=] {
                    for (size_t i = 0; i < 100; i++) {
                        expect(that % fabs(fsine_results[i] - sin_results[i]) < 1e-6);
                        expect(that % fabs(dsine_results[i] - sin_results[i]) < 1e-12);
                        expect(that % fabs(dsine_vec_results[i] - sin_results[i]) < 1e-12);
                        expect(that % fabs(dsine_unrolled_results[i] - sin_results[i]) < 1e-12);
                    }
                };
            };
        };
    };
};