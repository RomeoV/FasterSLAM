#include "add_observation_noise.h"  // import file to test
#include <cmath>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

    //! This seed was set in the add_observation_noise() function implemented in the yglee code,
    //! in order to be able to take the respective inputs and outputs of the function in that
    //! code, and test if we can exactly reproduce the results in our code.
    srand(1994);

    "add_observation_noise"_test = [] {
        given("I have some arguments") = [] {
            Vector2d z[2] = {{25.77444969441645, -1.4733774573801159}, {25.276107769232738, 0.13836675004531551}};
            const int zlen = 2;
            double R[4] = {0.010000000000000002, 0, 0, 0.00030461741909055634};
            const int addnoise = 1;

            when("I call add_observation_noise(z, zlen, R, addnoise)") = [&] {
                add_observation_noise(z, zlen, R, addnoise);

                then("This is equal with the actual result") = [&] {
                    Vector2d actual_z[2] = {{25.791968854179345, -1.4749004792119464},
                                            {25.206177150152218, 0.14779058881507551}};
                    for (int i = 0; i < zlen; i++) {    
                        expect(fabs(z[i][0] - actual_z[i][0]) < 1e-10);
                        expect(fabs(z[i][1] - actual_z[i][1]) < 1e-10);
                    }
                };
            };
        };
    };

}
