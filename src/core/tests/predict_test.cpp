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
    /*
    "predict"_test = [] {
        given("I have a particle at position (x=1,y=1,theta=pi/2)") = [] {
            auto is_close = [](double lhs, double rhs) {
                return  fabs(lhs - rhs) < 1e-4;
            };

            auto is_approx = [](double lhs, double rhs) {
                return  fabs(lhs - rhs) < 2*1e-1;
            };

            SWITCH_PREDICT_NOISE = 0;
            Vector3d initial_pos = {1, 1, M_PI/2};
            Particle p;
            initParticle(&p, 10, initial_pos);

            when("I give a constant control (steering angle) and calculate the time it takes to drive half a circle") = [&] {
                double Q[4] = {0, 0, 0, 0};

                const double V = 1; //random()*1.0/RAND_MAX;
                const double G = M_PI/4.;

                auto time_per_half_circle = [V] (double G) {  // note that static variables don't need to be captured
                    return 1 * (M_PI * WHEELBASE)/(V * sin(G));
                };

                double T = time_per_half_circle(G);
                int k_per_half_circle = int( T * V / 0.01);
                const double dt = T / k_per_half_circle;  // make sure T is divisable by dt, and one timestep make about 1cm progress

                then("Before moving, it's in the starting position") = [&] {
                    expect(that % is_close(p.xv[0],  initial_pos[0]) == true) << "p.xv[0] = " << p.xv[0];
                    expect(that % is_close(p.xv[1],  initial_pos[1]) == true) << "p.xv[1] = " << p.xv[1];
                    expect(that % is_close(p.xv[2], initial_pos[2]) == true) << "p.xv[2] = " << p.xv[2];
                };

                for (size_t i = 0; i < k_per_half_circle; i++) {
                    predict(&p, V, G, Q, WHEELBASE, dt);
                }
                then("After half a circle, the angle is opposite and the y values is the same again (but not he x value)!") = [&] {
                    expect(that % is_close(p.xv[0],  initial_pos[0]) == false) << "p.xv[0] = " << p.xv[0];
                    expect(that % is_approx(p.xv[1],  initial_pos[1]) == true ) << "p.xv[1] = " << p.xv[1];
                    expect(that % is_close(p.xv[2], -initial_pos[2]) == true ) << "p.xv[2] = " << p.xv[2];
                };

                for (size_t i = 0; i < k_per_half_circle; i++) {
                    predict(&p, V, G, Q, WHEELBASE, dt);
                }
                then("After another half a circle, x, y and the angle is the same as initial position") = [&] {
                    expect(that % is_close(p.xv[0], initial_pos[0]) == true);
                    expect(that % is_close(p.xv[1], initial_pos[1]) == true) << "p.xv[1] = " << p.xv[1];
                    expect(that % is_close(p.xv[2], initial_pos[2]) == true) << "p.xv[2] = " << p.xv[2];
                };

                for (size_t i = 0; i < k_per_half_circle; i++) {
                    predict(&p, V, G, Q, WHEELBASE, dt);
                }
                then("After one and a half circle, it's in the same pos as during the first test") = [&] {
                    expect(that % is_close(p.xv[0],  initial_pos[0]) == false) << "p.xv[0] = " << p.xv[0];
                    expect(that % is_approx(p.xv[1], initial_pos[1]) == true) << "p.xv[1] = " << p.xv[1];
                    expect(that % is_close(p.xv[2], -initial_pos[2]) == true) << "p.xv[2] = " << p.xv[2];
                };
            };

            delParticleMembers(&p);
        };
    };
    */
}
