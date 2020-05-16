#include "KF_cholesky_update.h"  // import file to test
#include <cmath>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

    "KF_cholesky_update"_test = [] {
        given("I have the arguments x, P, v, R, H") = [] {

            double x[2] = {3.2403905331533212, -25.689432087069857};
            double P[4] = {0.20063369668655512, 0.018909593226709744, 0.018909593226709744, 0.011875705723671498};
            double v[2] = {-0.017001037783700212, -0.010645013219889199};
            double R[4] = {0.010000000000000002, 0, 0, 0.00030461741909055634};
            double H[4] = {0.073819203427568675, -0.99727164063023443, 0.038893076096335785, 0.0028789105989867826};
            
            when("I call KF_cholesky_update(x, P, v, R, H)") = [&] {
                
                KF_cholesky_update(x, P, v, R, H);

                then("I get the updated values of x and P I want") = [=] {
                    double actual_x[2] = {3.1065907987134258, -25.693760147763445};
                    double actual_P[4] = {0.098876426893456063, 0.0071308261313384278, 0.0071308261313384278, 0.0055235811468849847};
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
};
