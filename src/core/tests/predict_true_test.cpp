#include "predict_true.h"  // import file to test
#include <math.h>

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    "predict_true"_test = [] {
        given("I have the arguments V, G, WB, dt, xv") = [] {

            double V = 1.0;
            double G = M_PI/3.0; //60°
            double WB = 0.1;
            double dt = 1.0;
            double xv[3] = {0.0,0.0,0.0};

            when("I call predict_true(V, G, WB, dt, xv)") = [&] {
                
                predict_true(V, G, WB, dt, &xv);

                then("I get the updated values of xv I want") = [=] {

                    double exact_xv[3] = {0.5,sqrt(3.0)/2.0, 2.377068730665};
                    for (int i = 0; i < 3; i++) {
                        double error = fabs(xv[i] - exact_xv[i]);
                        expect(that % error < 1e-12);
                    }
                };
            };
        };
    };
};
