#include "linalg.h"  // import file to test
#include <cmath>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

  "Id_times_Id"_test = [] {
    given("I have two identity matrices") = [] {
      const size_t N = 2;
      double A[N*N] = {1., 0., 0., 1.};
      double B[N*N] = {1., 0., 0., 1.};

      when("I multiply them") = [=] {
        double C[N*N];
        mul(A, B, N, N, N, C);
        then("They are still the identity matrix") = [C, N] {
            for (size_t r = 0; r < N; r++) {
              for (size_t c = 0; c < N; c++) {
                expect(fabs(C[r*N+c] - (r == c)) < 1e-14);
              }
            }
        };
      };
    };
  };


  "M_times_M_inv"_test = [] {
    given("I have a random matrix M") = [] {
      const size_t N = 2;
      double A[N*N];
      fill_rand(A, 4, -1, 1);

      when("I take the inverse M_inv") = [=] {
        double B[N*N];
        inv_2x2(A, B);

        then("M*M_inv and M_inv*M are the identity matrix") = [=] {
          double C_lhs[N*N];
          double C_rhs[N*N];
          mul(A, B, N, N, N, C_rhs);
          mul(B, A, N, N, N, C_lhs);

          "Compare M*M_inv to Id"_test = [C_rhs, N]() {
            for (size_t r = 0; r < N; r++) {
              for (size_t c = 0; c < N; c++) {
                expect(fabs(C_rhs[r*N+c] - (r == c)) < 1e-14) << "r[" << r << "], c[" << c << "]: Value is " << C_rhs[r*N+c];
              }
            }
          };

          "Compare M_inv*M to Id"_test = [C_lhs, N]() {
            for (size_t r = 0; r < N; r++) {
              for (size_t c = 0; c < N; c++) {
                expect(fabs(C_lhs[r*N+c] - (r == c)) < 1e-14) << "r[" << r << "], c[" << c << "]: Value is " << C_lhs[r*N+c];
              }
            }
          };
        };
      };
    };
  };

};

