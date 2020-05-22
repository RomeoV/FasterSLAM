#include "particle.h"

/*!
    Predict step of the algorithm, predicts the state of a single particle.
    Uses multivariate_gauss.
    @param[out] particle Particle whose position is predicted.
    @param[in]  V        velocity that has been generated by add_control_noise
    @param[in]  G        steer angle in (-pi,pi] that has been modified by add_control_noise
    @param[in]  Q        matrix of control noises as in configfile.cpp
    [REMOVED] param[in]  WB        wheelbase, float, metres. In configfile.cpp
    @param[in]  dt        seconds, time interval between control signals, float. In configfile.cpp
    @param[in]  addrandom        0/1, if sampling from predict noise, for fs1 usually true/1
 */
void predict(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt);
void predict_base(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt);
void predict_active(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt);

double predict_base_flops(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt);
double predict_base_memory(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt);