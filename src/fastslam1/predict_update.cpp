#include "predict_update.h"

#include "particle.h"

#include "compute_steering.h"
#include "predict_true.h"
#include "add_control_noise.h"
#include "predict.h"
#include "configfile.h"
void predict_update(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    predict_update_base(wp, N_waypoints, V, Q, dt, xtrue, iwp, G, particles);                   
}

void predict_update_base(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    compute_steering_base(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true_base(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2];
    add_control_noise_base(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step	
    for (size_t i = 0; i < NPARTICLES; i++) {
        predict_base(&particles[i], VnGn[0], VnGn[1], Q, WHEELBASE, dt);
    }
}

void predict_update_active(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    compute_steering(xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, iwp, G);
            if ( *iwp == -1 && NUMBER_LOOPS > 1 ) {
                *iwp = 0;
                NUMBER_LOOPS--;
            }
    predict_true(V, *G, WHEELBASE, dt, xtrue);

    // add process noise
    double VnGn[2];
    add_control_noise(V, *G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
    // Predict step	
    for (size_t i = 0; i < NPARTICLES; i++) {
        predict(&particles[i], VnGn[0], VnGn[1], Q, WHEELBASE, dt);
    }
}