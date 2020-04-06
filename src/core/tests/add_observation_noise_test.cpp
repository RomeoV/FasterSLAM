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

    "add_observation_noise"_test = [] {
        given("I have some arguments") = [] {
            double *z;
            const int zlen;
            cMatrix2d R;
            const int addnoise = 1;

            when("I call add_observation_noise(z, zlen, R, addnoise)") = [&] {
                add_observation_noise(z, zlen, R, addnoise);

                then("This is equal with the actual result") = [&] {
                    double *actual_z;
                    for (int i = 0; i < zlen; i++) {    
                        expect(fabs(z[i*2+0] - actual_z[i*2+0]) < 1e-10);
                        expect(fabs(z[i*2+1] - actual_z[i*2+1]) < 1e-10);
                    }
                };
            };
        };
    };

};
