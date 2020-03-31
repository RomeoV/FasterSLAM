#include "multivariate_gauss.h"  // import file to test
#include <cmath>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

    //! This seed was set in the multivariate_gauss() function implemented in the yglee code,
    //! in order to be able to take the respective inputs and outputs of the function in that
    //! code, and test if we can exactly reproduce the results in our code.
    srand(1994);

    "multivariate_gauss"_test = [] {
        given("I have the arguments x and P") = [] {
            double x[2] = {3.0, -0.008726646250000001};
            double P[4] = {0.089999999999999997, 0, 0, 0.0027415567718150069};

            when("I call multivariate_gauss(x, P, res)") = [=] {
                double res[2];
                multivariate_gauss(x, P, res);

                then("This is equal with the actual result") = [=] {
                    double actual_res[2] = {3.0525574792886885, -0.045342232735551199};
                    expect(fabs(res[0] - actual_res[0]) < 1e-10);
                    expect(fabs(res[1] - actual_res[1]) < 1e-10);
                };
            };
        };
    };

};
