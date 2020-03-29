#include "linalg.h"  // import file to test
#include <cmath>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

//!
//! Matrix multiplication test
//!
"Id_times_Id"_test = [] {
    given("I have two identity matrices") = [] {
        const size_t n = 2;
        double A[n*n] = {1., 0., 0., 1.};
        double B[n*n] = {1., 0., 0., 1.};

        when("I multiply them") = [=] {
            double C[n*n];
            mul(A, B, n, n, n, C);
            then("They are still the identity matrix") = [C, n] {
                for (size_t i = 0; i < n; i++) {
                    for (size_t j = 0; j < n; j++) {
                        expect(fabs(C[i*n+j] - (i == j)) < 1e-14);
                    }
                }
            };
        };
    };
};

//!
//! Matrix multiplication and 2x2 matrix inversion test
//!
"A_times_A_inv"_test = [] {
    given("I have a random matrix A") = [] {
        const size_t n = 2;
        double A[n*n];
        fill_rand(A, 4, -1, 1);

        when("I take the inverse Ainv") = [=] {
            double Ainv[n*n];
            inv_2x2(A, Ainv);

            then("Ainv*A and A*Ainv are the identity matrix") = [=] {
                double Il[n*n];
                double Ir[n*n];
                mul(Ainv, A, n, n, n, Il);
                mul(A, Ainv, n, n, n, Ir);

                "Compare Ainv*A to Id"_test = [Il, n]() {
                    for (size_t i = 0; i < n; i++) {
                        for (size_t j = 0; j < n; j++) {
                            expect(fabs(Il[i*n+j] - (i == j)) < 1e-14) << "i[" << i << "], j[" << j << "]: Value is " << Il[i*n+j];
                        }
                    }
                };

                "Compare A*Ainv to Id"_test = [Ir, n]() {
                    for (size_t i = 0; i < n; i++) {
                        for (size_t j = 0; j < n; j++) {
                            expect(fabs(Ir[i*n+j] - (i == j)) < 1e-14) << "i[" << i << "], j[" << j << "]: Value is " << Ir[i*n+j];
                        }
                    }
                };
            };
        };
    };
};

//!
//! Cholesky decomposition (2x2 matrix), matrix transposition and matrix multiplication test
//!
"L_times_Lt"_test = [] {
    given("I have an SPD matrix A") = [] {
        const size_t n = 2;
        double A[n*n] = {2., -1., -1., 2.};

        when("I compute the cholesky decomposition A = L*L^T") = [=] {
            double L[n*n], Lt[n*n];;
            llt_2x2(A, L);
            transpose(L, n, n, Lt);

            then("L*L^T is equal to A") = [=] {
                double LLt[n*n];
                mul(L, Lt, n, n, n, LLt);

                "Compare L*L^T to A"_test = [LLt, A, n]() {
                    for (size_t i = 0; i < n; i++) {
                        for (size_t j = 0; j < n; j++) {
                            expect(fabs(LLt[i*n+j] - A[i*n+j]) < 1e-14) << "i[" << i << "], j[" << j << "]: Value is " << LLt[i*n+j];
                        }
                    }
                };
            };
        };
    };
};

}

