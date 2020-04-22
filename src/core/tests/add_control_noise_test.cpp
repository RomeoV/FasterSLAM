#include "add_control_noise.h"  // import file to test
#include <cmath>
#include <sstream>
#include <algorithm>

#include "typedefs.h"
#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

"add_noise"_test = [] {
    // since this is random, it is a bit hard to test
    // for now test that it changed
    given("I have speed, a steering angle and a covariance matrix") = [] {
        double V = 10;
        double G = 0.5;
        double Q[4] = {0.1, 0.1, 0.3, 0.5};
        int addnoise = 1; // adding noise
        when("I add control noise") = [&] {
            double VnGn[2];
            add_control_noise(V, G, Q, addnoise, VnGn);
            then("I get the other speed and velocity than I had before") = [&] {
                bool different = ((VnGn[0] != V) && (VnGn[1] != G));
                expect(that % different == true);
            };
        };
    };
};

//!
//! Matrix multiplication test
//!
"don't add noise"_test = [] {
    // the function shouldn't do anything in this case
    given("I have speed, a steering angle and a covariance matrix") = [] {
        double V = 10;
        double G = 0.5;
        double Q[4] = {0.1, 0.1, 0.3, 0.5};
        int addnoise = 0; // not adding noise
        when("I do not add control noise") = [&] {
            double VnGn[2] = {0., 0.};
            add_control_noise(V, G, Q, addnoise, VnGn);
            then("I get the same speed and velocity I had before") = [&] {
                expect(that % VnGn[0] == 0.);
                expect(that % VnGn[1] == 0.);
            };
        };
    };
};

}

