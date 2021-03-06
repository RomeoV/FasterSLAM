#include "transform_to_global.h"  // import file to test
#include <cmath>
#include <sstream>
#include <algorithm>

#include "typedefs.h"
#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

"points"_test = [] {
    given("I have 2 points in 2D and a vector with a x/y/rot transform") = [] {
        Matrix2d p = {1., 1., 3., 1.};
        int size = 4;
        Vector3d b = {1., 1., M_PI*2};
        when("I transform the points") = [&] {
            transform_to_global(p, size, b);
            then("I get the same points moved 1 in and 1 in y direction") = [&] {
                expect(that % p[0] == 2);
                expect(that % p[1] == 2);
                expect(that % p[2] == 4);
                expect(that % p[3] == 2);
            };
        };
    };
};

}

