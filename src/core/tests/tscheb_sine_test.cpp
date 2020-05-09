#include "tscheb_sine.h"  // import file to test
#include "pi_to_pi.h"
#include <cmath>
#include <memory>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    "tscheb_sine"_test = [] {
        given("I have 99 random angles in [-10PI, 10PI]") = [] {
            double* angles = (double*)aligned_alloc(32, 100*sizeof(double));
            auto rand_angle = []{ return rand()*1.0/RAND_MAX * 4. * M_PI - 2. * M_PI; };
            for (size_t i = 0; i < 99; i++) { angles[i] = rand_angle(); }

            when("I calculate the cmath::sin, the tscheb_fsine and the tscheb_dsine") = [&] {
                double sin_results[99];
                double fsine_results[99];
                double dsine_results[99];
                double dsine_vec_results[99];
                double dsine_unrolled_results[99];
                double* dsine_avx_results = (double*)aligned_alloc(32, 100*sizeof(double));
                double* dsine_avx_unrolled_results = (double*)aligned_alloc(32, 100*sizeof(double));

                for (size_t i = 0; i < 99; i++) {
                    sin_results[i] = sin(angles[i]);
                    fsine_results[i] = tscheb_fsine(angles[i], false);
                    dsine_results[i] = tscheb_dsine(angles[i], false);
                }

                // normalize angles for methods that can only deal with normalized angles
                for (size_t i = 0; i < 99; i++) {
                    angles[i] = pi_to_pi(angles[i]);
                }
                tscheb_dsines(angles, 99, dsine_vec_results);
                tscheb_dsines_unrolled(angles, 99, dsine_unrolled_results);
                "aligned_avx_containers"_test = [&]{
                    expect(((unsigned long)angles & 15) == 0);
                    expect(((unsigned long)dsine_avx_results & 15) == 0);
                    expect(((unsigned long)dsine_avx_unrolled_results & 15) == 0);
                };
                tscheb_dsines_avx(angles, 99, dsine_avx_results);
                tscheb_dsines_avx_unrolled(angles, 99, dsine_avx_unrolled_results);

                then("I expect each approximation to be close to cmath::sin, up to some tolerance") = [=] {
                    for (size_t i = 0; i < 99; i++) {
                        expect(that % fabs(fsine_results[i] - sin_results[i]) < 1e-6);
                        expect(that % fabs(dsine_results[i] - sin_results[i]) < 1e-6);
                        expect(that % fabs(dsine_vec_results[i] - sin_results[i]) < 1e-6);
                        expect(that % fabs(dsine_unrolled_results[i] - sin_results[i]) < 1e-6);
                        expect(that % fabs(dsine_avx_results[i] - sin_results[i]) < 1e-6);
                        expect(that % fabs(dsine_avx_unrolled_results[i] - sin_results[i]) < 1e-6);
                    }
                };

                free(dsine_avx_results);
                free(dsine_avx_unrolled_results);
            };
            free(angles);
        };
    };
};