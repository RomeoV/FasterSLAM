#include "predict.h"  // import file to test
#include <math.h>
#include <random>
#include "configfile.h"

#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

/**
 * I went through the math and came up with a way to test.
 * The motion model is defined as follows:
 * \f[
 * \vec{x}(t) = [x(t), y(t), \alpha(t)]^T \\
 * 
 * \vec{x}'(t) = [V \cdot \cos(\alpha),
 *                V \cdot \sin(\alpha),
 *                V \cdot \frac{\sin(G)}{WB}]^T
 * \f]
 * 
 * Note in particular that \f$x'[2](t)\f$ is only dependent on the input \f$G(t)\f$.
 * This means with a given \f$G(t)\f$ we can analytically derive the solution of \f$x[2](t)\f$ and then derive either \f$x(t)\f$ or \f$y(t)\f$.
 * Let \f$G(t) = G\f$ be constant.
 * Then we can solve
 * \f[
 * \begin{align}
 * \alpha(t) &= \int_0^t V \cdot \sin(G) / WB \cdot dt \\
 *           &= t \cdot V \cdot \sin(G) / WB
 * \end{align}
 * \f]
 * 
 * Now we can derive \f$x(t)\f$:
 * \f[
 * \begin{align}
 * x(t) &= \int_0^t V \cdot cos(\alpha(t)) \cdot dt \\
 *      &= \int_0^t V \cdot cos(t \cdot V \cdot \sin(G) / WB) \cdot dt \\
 *      & \text{let } C = V \cdot \sin(G) / WB \text{ and substitute } s = t \cdot C \\
 *      &= \int_0^{t \cdot C} cos(s) \cdot ds \cdot 1/C \\
 *      &= \sin(s)[t \cdot C; 0] \cdot 1/C \\
 *      &= \sin(t \cdot C) \cdot 1/C
 * \end{align}
 * \f]
 * Note that with a constant G, we expect the car to drive in a circle.
 * Thus, assuming the car starts at \f$[0,0,0]^T\f$, let's calculate the time when x passes through zero again, i.e. either a half or a full circle has been made:
 * \f[
 * \begin{align}
 * & x(t) \stackrel{!}{=} 0 \\
 * \Leftrightarrow& \sin(t \cdot C) \cdot 1/C = 0 \\
 * \Leftrightarrow& t \cdot C = k \cdot pi \\
 * \Leftrightarrow& t = k \cdot \frac{\pi \cdot (WB)}{V \cdot \sin(G)} \text{ with } k \in\{1,2,\dots\}
 * \end{align}
 * \f]
 * Note that the derivation for \f$y\f$ is basically the same.
 * Thus, we can test the following:
 * Speficy some \f$V, WB, G\f$, calculate the time t for a half, full and double cirlce and each time check \f$x\f$ and \f$y\f$ if they are close to the expected value.
 * Note that because some basic Euler scheme is used, the precision will likely not be super high...
 */

int main() {
    "predict"_test = [] {
        given("I have a particle at position (x=1,y=1,theta=pi/2)") = [] {
            SWITCH_PREDICT_NOISE = 0;
            double initial_pos[3] = {1, 1, M_PI/2};
            Particle p;
            p.xv[0] = 1; 
            p.xv[1]=1;
            p.xv[2]= M_PI/2;

            when("I give control inputs for a rectangle") = [=]() mutable {
                double Q[4] = {0, 0, 0, 0};

                const double V = random()*1.0/RAND_MAX;
                const double dt = random()*1.0/RAND_MAX;
                double intermediate_pos[3] = {1-2*V*dt, 1 - V*dt, -M_PI/2};

                //predict(&p, 2*V, M_PI/2 + 0*M_PI, Q, dt);
                //predict(&p, 1*V, M_PI/2 + 2*M_PI, Q, dt);  // +2PI to check invariance to values outside of (-PI,PI]
                then("After half way, I'm in a different place") = [&] {
                    "position"_test = [&](size_t i) {
                        expect(that % fabs(intermediate_pos[i] - p.xv[i]) < 1e-14) << p.xv[i] << "!= " << intermediate_pos[i];
                    } | std::vector{0,1,2};
                };

                //predict(&p, 2*V, M_PI/2 + 4*M_PI, Q, dt);
                //predict(&p, 1*V, M_PI/2 - 2*M_PI, Q, dt);
                then("In the end I end up in the same place as before") = [&] {
                    "position"_test = [&](size_t i) {
                        expect(that % fabs(initial_pos[i] - p.xv[i]) < 1e-14) << p.xv[i];
                    } | std::vector{0,1,2};
                };
            };

            when("I give control inputs for a pentagon") = [=]() mutable {
                double Q[4] = {0, 0, 0, 0};

                const double V = random()*1.0/RAND_MAX;
                const double dt = random()*1.0/RAND_MAX;
                //predict(&p, 1*V, M_PI*2./5 + 0*M_PI, Q, dt);
                //predict(&p, 1*V, M_PI*2./5 + 6*M_PI, Q, dt);
                //predict(&p, 1*V, M_PI*2./5 + 2*M_PI, Q, dt);
                then("After 3/5 of the way, I'm in a different place") = [&] {
                    "position"_test = [&](size_t i) {
                        expect(that % fabs(initial_pos[i] - p.xv[i]) > 1e-14) << p.xv[i];
                    } | std::vector{0,1,2};
                };

                //predict(&p, 1*V, M_PI*2./5 - 8*M_PI, Q, dt);
                //predict(&p, 1*V, M_PI*2./5 - 2*M_PI, Q, dt);
                then("I end up in the same place as before") = [&] {
                    "position"_test = [&](size_t i) {
                        expect(that % fabs(initial_pos[i] - p.xv[i]) < 1e-14) << p.xv[i];
                    } | std::vector{0,1,2};
                };
            };
        };
    };
}
