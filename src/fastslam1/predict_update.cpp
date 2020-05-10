#include "predict_update.h"

#include "particle.h"

#include "compute_steering.h"
#include "predict_true.h"
#include "add_control_noise.h"
#include "predict.h"
#include "configfile.h"
void predict_update(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
    predict_update_base(wp, N_waypoints, V, Q, dt, N, xtrue, iwp, G, particles);                   
}


/*****************************************************************************
 * PERFORMANCE STATUS (N=NPARTICLES)
 * Work, best: 20 (pred_true) + 22(com_steer.) + N*20 (predict) = 42 + N*20 flops
 * Work, worst: 29 (pred_true) +17 (con_noise) + 37(com_steer.) + N*46 (predict) = 83 + N*46 flops
 * 
 * 
 * #work best, det: 1 atan2 + 2 pow_2 + 3*N+3 sin/cos + 9*N+9+3 mults + 3*N+3+8 adds + 1*N+1+1 neg + 4*N+4+7 fl-comp + 1*N+1 div  
 * #work, worst, det: 1 atan2 + 2 pow_2 + 3*N+3+2 neg + 17*N+12+9 mults + 4*N+2+1 div + 1*N+1+1 floor + 4*N+4+11 fl-comp + 3*N+3 sin + 12*N+5+10 adds + 2*N sqrt 
 * #Work best, detailed: 
 * Memory moved: TBD
 * Cycles: 48000 with N = 100
 * Performance: 0.04
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void predict_update_base(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
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
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles) {
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
    for (size_t i = 0; i < N; i++) {
        predict(&particles[i], VnGn[0], VnGn[1], Q, WHEELBASE, dt);
    }
}