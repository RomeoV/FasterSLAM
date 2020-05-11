#include "linalg.h"  // import file to test
#include <cmath>
#include <sstream>
#include <algorithm>

#include "typedefs.h"
#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

"v*A*vT"_test = [] {
    given("I have a 3-dim vector v of ones and a 3x3 matrix M with all ones") = [] {
        Vector3d v = {1., 1., 1.};
        Matrix3d M = {1., 1., 1.,
                      1., 1., 1.,
                      1., 1., 1.};
        when("I multiply v*M*v.T") = [&] {
            Vector3d M_vT;
            double v_M_vT;
            mul(M, v, 3, 3, 1, M_vT);
            mul(v, M_vT, 1, 3, 1, &v_M_vT);
            then("I get a double with value 9") = [&] {
                expect(that % fabs(v_M_vT - 9) < 1e-14);
            };
        };
    };
};

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

                "Compare Ainv*A to Id"_test = [Il, n] {
                    for (size_t i = 0; i < n; i++) {
                        for (size_t j = 0; j < n; j++) {
                            expect(fabs(Il[i*n+j] - (i == j)) < 1e-14) << "i[" << i << "], j[" << j << "]: Value is " << Il[i*n+j];
                        }
                    }
                };

                "Compare A*Ainv to Id"_test = [Ir, n] {
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

                "Compare L*L^T to A"_test = [LLt, A, n] {
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

//! Dummy fill test
"fill"_test = [] {
    given("I have a static array initialized to zero and a variable val set to some value") = [] {
        const int n = 4;
        double x[n] = { };
        const double val = 2.0;

        when("I call fill(x, n, val)") = [&] {
            fill(x, n, val);

            then("All elements of x are set to the value of val") = [&] {
                "Compare x[i]s to val"_test = [&] {
                    for (int i = 0; i < n; i++) {
                        expect( x[i] == val ) << "x[" << i << "] = " << x[i] << " != val = " << val;
                    }
                };
            };
        };
    };
};

"print"_test = [] {
    given("I have a 3x2 matrix with running numbers 1..6") = [] {
        const size_t n = 3, m = 2;
        double x[n*m] = {1, 2, 3, 4, 5, 6};

        when("I print it") = [=] {
            std::stringstream ss;
            print(x, n, m, ss);
            std::string str = ss.str();

            then("I find all numbers after each other, seperated by commas, whitespace or line breaks") = [=] {
                auto is_whitespace_or_comma = [](char c){return c == ' ' or c == '\t' or c == ',';};

                auto iter = str.begin();
                for (size_t r = 0; r < n; r++) {
                    for (size_t c = 0; c < m-1; c++) {
                        iter = std::find(iter, str.end(), '0'+(r*m)+c+1);
                        expect(iter != str.end()) << "Couldn't find " << (r*m)+c+1;

                        iter = std::find_if(iter, str.end(), is_whitespace_or_comma);
                        expect(iter != str.end()) << "Couldn't find whitespace";
                    }

                    iter = std::find(iter, str.end(), '0'+(r+1)*m);
                    expect(iter != str.end());

                    iter = std::find(iter, str.end(), '\n');
                    expect(iter != str.end()) << "Couldn't find newline";
                }
            };
        };
    };
};

"transpose_2x2"_test = [] {
    given("I have a matrix (row-major)") = [] {
        double A[4] = {1., 2., 3., 4.};
        when("I transpose it") = [=] {
            double T[4];
            transpose_2x2(A, T);
            then("I get the transpose (col-major)") = [=] {
                double T_actual[4] = {1., 3., 2., 4.};
                for (size_t i = 0; i < 4; i++) {
                    expect(fabs(T[i] - T_actual[i]) < 1e-16);
                }
            };
        };
    };
};

"stranspose_2x2"_test = [] {
    given("I have a matrix (row-major)") = [] {
        double A[4] = {1., 2., 3., 4.};
        when("I transpose it") = [&] {
            stranspose_2x2(A);
            then("I get the transpose (col-major)") = [=] {
                double T_actual[4] = {1., 3., 2., 4.};
                for (size_t i = 0; i < 4; i++) {
                    expect(fabs(A[i] - T_actual[i]) < 1e-16);
                }
            };
        };
    };
};


//!
//! Matrix - Matrix multiplication test (2x2)
//!
"mm_2x2"_test = [] {
    given("I have two identity matrices") = [] {
        double A[4] = {1., 0., 0., 1.};
        double B[4] = {1., 0., 0., 1.};

        when("I multiply them") = [=] {
            double C[4];
            mm_2x2(A, B, C);
            then("They are still the identity matrix") = [=] {
                for (size_t i = 0; i < 2; i++) {
                    for (size_t j = 0; j < 2; j++) {
                        expect(fabs(C[i*2+j] - (i == j)) < 1e-16);
                    }
                }
            };
        };
    };
    given("I have two matrices") = [] {
        double A[4] = {1., 2., 3., 4.};
        double B[4] = {5., 6., 7., 8.};

        when("I multiply them") = [=] {
            double C[4], D[4];
            mm_2x2(A, B, C);
            mul(A, B, 2, 2, 2, D);
            then("I get the correct result ( according to mul() )") = [=] {
                for (size_t i = 0; i < 2; i++) {
                    for (size_t j = 0; j < 2; j++) {
                        expect(fabs(C[i*2+j] - D[i*2+j]) < 1e-16);
                    }
                }
            };
        };
    };
    given("I have two matrices") = [] {
        double A[4] = {1., 2., 3., 4.};
        double B[4] = {5., 6., 7., 8.};

        when("I multiply A * B^T") = [&] {
            double C[4], D[4];
            mmT_2x2(A, B, C);
            stranspose_2x2(B);
            mul(A, B, 2, 2, 2, D);
            then("I get the correct result ( according to stranspose() and mul() )") = [=] {
                for (size_t i = 0; i < 2; i++) {
                    for (size_t j = 0; j < 2; j++) {
                        expect(fabs(C[i*2+j] - D[i*2+j]) < 1e-16);
                    }
                }
            };
        };
    };
    given("I have three matrices") = [] {
        double A[4] = {1., -5., 3.14, 0.};
        double B[4] = {12., 45., 56., 16.};
        double C[4] = {133., 0.2, 1.2, 5.6}; // will be modified inside mmadd_2x2
        when("I multiply-add them") = [&] {
            double D[4];
            mul(A, B, 2, 2, 2, D); // D = A*B
            add(D, C, 4, D);       // D = D + C
            mmadd_2x2(A, B, C);    // C = A*B + C 
            then("I get the correct matrix ( according to mul() and add() )") = [=] {
                for (size_t i = 0; i < 4; i++) {
                    expect(fabs(C[i] - D[i]) < 1e-16);
                }
            };
        };
    };
};

//!
//! Matrix - Vector multiplication test (2x2)
//!
"mv_2x2"_test = [] {
    given("I have an identity matrix and a vector") = [] {
        double A[4] = {1., 0., 0., 1.};
        double b[2] = {12., 45.};

        when("I multiply them") = [=] {
            double c[2];
            mv_2x2(A, b, c);
            then("I get the exact same vector") = [=] {
                for (size_t i = 0; i < 2; i++) {
                    expect(fabs(c[i] - b[i]) < 1e-16);
                }
            };
        };
    };
    given("I have a matrix and a vector") = [] {
        double A[4] = {1., -5., 3.14, 0.};
        double b[2] = {12., 45.};

        when("I multiply them") = [=] {
            double c[2], d[2];
            mv_2x2(A, b, c);
            mul(A, b, 2, 2, 1, d);
            then("I get the correct vector ( according to mul() )") = [=] {
                for (size_t i = 0; i < 2; i++) {
                    expect(fabs(c[i] - d[i]) < 1e-16);
                }
            };
        };
    };
    given("I have a matrix and two vectors") = [] {
        double A[4] = {1., -5., 3.14, 0.};
        double b[2] = {12., 45.};
        double c[2] = {133., 0.2}; // will be modified inside mvadd_2x2
        when("I multiply-add them") = [&] {
            double d[2];
            mul(A, b, 2, 2, 1, d); // d = A*b
            add(d, c, 2, d);       // d = d + c
            mvadd_2x2(A, b, c);    // c = A*b + c 
            then("I get the correct vector ( according to mul() and add() )") = [=] {
                for (size_t i = 0; i < 2; i++) {
                    expect(fabs(c[i] - d[i]) < 1e-16);
                }
            };
        };
    };
};

}

