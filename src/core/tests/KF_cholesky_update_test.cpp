#include "KF_cholesky_update.h"  // import file to test
#include <cmath>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
#ifdef KF_YGLEE
    "KF_cholesky_update"_test = [] {
        given("I have the arguments x, P, v, R, H") = [] {

            double x[2] __attribute__((aligned(32))) = {3.2403905331533212, -25.689432087069857};
            //double x[2] = {3.227460886446243, -25.613382543676146};
            double P[4] __attribute__((aligned(32))) = {0.20063369668655512, 0.018909593226709744, 0.018909593226709744, 0.011875705723671498};
            //double P[4] = {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843};
            double v[2] __attribute__((aligned(32))) = {-0.017001037783700212, -0.010645013219889199};
            //double v[2] = {0.128762949830296, 0.019814250533567};
            double R[4] __attribute__((aligned(32))) = {0.010000000000000002, 0, 0, 0.00030461741909055634};
            //double R[4] = {0.010000000000000, 0.0, 0.0, 0.000304617419787};
            double H[4] __attribute__((aligned(32))) = {0.073819203427568675, -0.99727164063023443, 0.038893076096335785, 0.0028789105989867826};
            //double H[4] = {0.075431770172036, -0.997150965525639, 0.038902059641499, 0.002942835461779};
            
            when("I call KF_cholesky_update(x, P, v, R, H)") = [&] {
                
                KF_cholesky_update(x, P, v, R, H);

                then("I get the updated values of x and P I want") = [=] {
                    double actual_x[2] = {3.1065907987134258, -25.693760147763445};
                    //double actual_x[2] = {3.470171202213126, -25.656742169761873};
                    double actual_P[4] = {0.098876426893456063, 0.0071308261313384278, 0.0071308261313384278, 0.0055235811468849847};
                    //double actual_P[4] = {0.099441292170602, 0.008341518741166, 0.008341518741166, 0.005737523050281};
                    for (int i = 0; i < 2; i++) {
                        expect(fabs(x[i] - actual_x[i]) < 1e-10) << x[i] << " != " << actual_x[i];
                    }
                    for (int i = 0; i < 4; i++) {
                        expect(fabs(P[i] - actual_P[i]) < 1e-10) << P[i] << " != " << actual_P[i];
                    }
                };
            };
        };
    };
#else
    "KF_cholesky_update"_test = [] {
        given("I have the arguments x, P, v, R, H") = [] {

            //double x[2] __attribute__((aligned(32))) = {3.2403905331533212, -25.689432087069857};
            double x[2] __attribute__((aligned(32))) = {3.227460886446243, -25.613382543676146};
            //double P[4] __attribute__((aligned(32))) = {0.20063369668655512, 0.018909593226709744, 0.018909593226709744, 0.011875705723671498};
            double P[4] __attribute__((aligned(32))) = {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843};
            //double v[2] __attribute__((aligned(32))) = {-0.017001037783700212, -0.010645013219889199};
            double v[2] __attribute__((aligned(32))) = {0.128762949830296, 0.019814250533567};
            //double R[4] __attribute__((aligned(32))) = {0.010000000000000002, 0, 0, 0.00030461741909055634};
            double R[4] __attribute__((aligned(32))) = {0.010000000000000, 0.0, 0.0, 0.000304617419787};
            //double H[4] __attribute__((aligned(32))) = {0.073819203427568675, -0.99727164063023443, 0.038893076096335785, 0.0028789105989867826};
            double H[4] __attribute__((aligned(32))) = {0.075431770172036, -0.997150965525639, 0.038902059641499, 0.002942835461779};
            
            when("I call KF_cholesky_update(x, P, v, R, H)") = [&] {
                
                KF_cholesky_update(x, P, v, R, H);

                then("I get the updated values of x and P I want") = [=] {
                    //double actual_x[2] = {3.1065907987134258, -25.693760147763445};
                    double actual_x[2] = {3.470171202213126, -25.656742169761873};
                    //double actual_P[4] = {0.098876426893456063, 0.0071308261313384278, 0.0071308261313384278, 0.0055235811468849847};
                    double actual_P[4] = {0.099441292170602, 0.008341518741166, 0.008341518741166, 0.005737523050281};
                    for (int i = 0; i < 2; i++) {
                        expect(fabs(x[i] - actual_x[i]) < 1e-10) << x[i] << " != " << actual_x[i];
                    }
                    for (int i = 0; i < 4; i++) {
                        expect(fabs(P[i] - actual_P[i]) < 1e-10) << P[i] << " != " << actual_P[i];
                    }
                };
            };
        };
    };
#ifdef __AVX2__
    "KF_cholesky_update"_test = [] {
        given("I have the arguments x, P, v, R, H") = [] {

            double x[4] __attribute__((aligned(32))) = {3.227460886446243, -25.613382543676146, 3.227460886446243, -25.613382543676146};
            double P[4] __attribute__((aligned(32))) = {0.199855490073439, 0.019180472296076, 0.019180472296076, 0.011937739684843};
            double v[4] __attribute__((aligned(32))) = {0.128762949830296, 0.019814250533567, 0.128762949830296, 0.019814250533567};
            double R[4] __attribute__((aligned(32))) = {0.010000000000000, 0.0, 0.0, 0.000304617419787};
            double H[4] __attribute__((aligned(32))) = {0.075431770172036, -0.997150965525639, 0.038902059641499, 0.002942835461779};
            
            __m256d x0x2 = _mm256_load_pd( x );
            __m256d x1x3 = _mm256_load_pd( x );
            __m256d P0   = _mm256_load_pd( P );
            __m256d P1   = _mm256_load_pd( P );
            __m256d P2   = _mm256_load_pd( P );
            __m256d P3   = _mm256_load_pd( P );
            __m256d v0v2 = _mm256_load_pd( v );
            __m256d v1v3 = _mm256_load_pd( v );
            __m256d RR   = _mm256_load_pd( R );
            __m256d H0   = _mm256_load_pd( H );
            __m256d H1   = _mm256_load_pd( H );
            __m256d H2   = _mm256_load_pd( H );
            __m256d H3   = _mm256_load_pd( H );
            when("I call KF_cholesky_update_unrolled4_avx()") = [&] {

                KF_cholesky_update_unrolled4_avx(&x0x2, &x1x3, &P0, &P1, &P2, &P3, v0v2, v1v3, RR, H0, H1, H2, H3);

                double x0x2_[4], x1x3_[4], P0_[4], P1_[4], P2_[4], P3_[4]; 
                _mm256_store_pd( x0x2_, x0x2 );
                _mm256_store_pd( x1x3_, x1x3 );
                _mm256_store_pd( P0_, P0 );
                _mm256_store_pd( P1_, P1 );
                _mm256_store_pd( P2_, P2 );
                _mm256_store_pd( P3_, P3 );

                then("I get the updated values of x and P I want") = [=] {
                    double actual_x[4] = {3.470171202213126, -25.656742169761873, 3.470171202213126, -25.656742169761873};
                    double actual_P[4] = {0.099441292170602, 0.008341518741166, 0.008341518741166, 0.005737523050281};
                    for (int i = 0; i < 4; i++) {
                        expect(fabs(x0x2_[i] - actual_x[i]) < 1e-10) << x0x2_[i] << " != " << actual_x[i];
                        expect(fabs(x1x3_[i] - actual_x[i]) < 1e-10) << x1x3_[i] << " != " << actual_x[i];
                    }
                    for (int i = 0; i < 4; i++) {
                        expect(fabs(P0[i] - actual_P[i]) < 1e-10) << P0[i] << " != " << actual_P[i];
                        expect(fabs(P1[i] - actual_P[i]) < 1e-10) << P1[i] << " != " << actual_P[i];
                        expect(fabs(P2[i] - actual_P[i]) < 1e-10) << P2[i] << " != " << actual_P[i];
                        expect(fabs(P3[i] - actual_P[i]) < 1e-10) << P3[i] << " != " << actual_P[i];
                    }
                };
            };
        };
    };
#endif
#endif
};
