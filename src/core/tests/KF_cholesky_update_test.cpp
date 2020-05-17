#include "KF_cholesky_update.h"  // import file to test
#include <cmath>
#include <vector>
#include <functional>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

    std::vector<std::pair<std::string, std::function<void (Vector2d, Matrix2d, cVector2d, cMatrix2d, cMatrix2d)>>>
    KF_functions = {
        {"KF_cholesky_update", KF_cholesky_update},
        {"KF_cholesky_update_base", KF_cholesky_update_base},
        //{KF_cholesky_update_v1},
        {"KF_cholesky_update_v2", KF_cholesky_update_v2}
    };
            
    "KF_cholesky_update"_test = [] (auto NamedFunc) {
        given("I have the arguments x, P, v, R, H") = [&] {
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

            when("I call " + NamedFunc.first + "(x, P, v, R, H)") = [&] {
                
                double x[2] __attribute__((aligned(32))) = {3.2403905331533212, -25.689432087069857};
                //double x[2] = {3.227460886446243, -25.613382543676146};
                NamedFunc.second(x, P, v, R, H);

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
    } | KF_functions;
};
