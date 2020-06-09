#include "trigonometry.h"  // import file to test

#include <vector>  // used for test input
#include "ut.hpp"  // import functionality and namespaces from single header file
using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`
#include <cmath>
#include <random>


int main() {
//#ifdef _  // this test fails...
    "atan2"_test = [] {
        given("A vector of 100 random x and y coords") = [] {
            std::cout << "atan2 test" << std::endl;
            const size_t N = 100;
            double xs[N], ys[N], cmath_results[N], approx_results[N];
            when("I calculate the cmath atan2 and the atan2-approximation") = [&] {
                for (size_t i = 0; i < N; i++) {
                    xs[i] = random() * 1.0/RAND_MAX * 10 - 5;
                    ys[i] = random() * 1.0/RAND_MAX * 10 - 5;
                    cmath_results[i] = atan2(ys[i], xs[i]);

                    double result = atan2_approximation1(ys[i], xs[i]);
                    approx_results[i] = result;
                }
                then("I expect the calculations to be rougly equal") = [&] {
                    for (size_t i = 0; i < N; i++) {
                        // only accurate to 1e-2*5
                        expect(that % fabs(approx_results[i] - cmath_results[i]) < 1e-2*5);
                    }
                };
            };
        };
    };

    "atan2-SIMD"_test = [] {
        given("A vector of 100 random x and y coords") = [] {
            std::cout << "atan2 test" << std::endl;
            const size_t N = 100;
            double xs[N], ys[N], cmath_results[N], approx_results[N];
            when("I calculate the cmath atan2 and the atan2-approximation") = [&] {
                for (size_t i = 0; i < N; i++) {
                    xs[i] = random() * 1.0/RAND_MAX * 10 - 5;
                    ys[i] = random() * 1.0/RAND_MAX * 10 - 5;

                    cmath_results[i] = atan2(ys[i], xs[i]);
                    double tmp[4] = {0.0,0.0,0.0,0.0};
                    __m256d ysimd = _mm256_set_pd(ys[i], ys[i], ys[i], ys[i]);
                    __m256d xsimd = _mm256_set_pd(xs[i], xs[i], xs[i], xs[i]);
                    __m256d result = atan2_approximation2(ysimd, xsimd);
                    _mm256_store_pd(tmp, result);
                    approx_results[i] = tmp[0];
                }
                then("I expect the calculations to be rougly equal") = [&] {
                    for (size_t i = 0; i < N; i++) {
                        // only accurate to 1e-2*5
                        expect(that % fabs(approx_results[i] - cmath_results[i]) < 1e-2*5);
                    }
                };
            };
        };
    };
//#endif
}
