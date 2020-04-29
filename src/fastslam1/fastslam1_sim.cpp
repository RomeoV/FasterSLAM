#include <iostream>
#include <math.h>
#include <string.h>
#include "add_control_noise.h"
#include "add_feature.h"
#include "add_observation_noise.h"
#include "compute_weight.h"
#include "data_associate_known.h"
#include "fastslam1_sim.h"
#include "feature_update.h"
#include "get_observations.h"
#include "linalg.h"
#include "predict.h"
#include "resample_particles.h"
#include "fastslam1_utils.h"
#include "predict_update.h"
#include "observe_update.h"

void fastslam1_sim( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_) 
{
    const size_t N_features = lm_rows;
    const size_t N_waypoints = wp_rows;

    Particle *particles;
    double *weights;
    Vector3d xtrue   = {0,0,0};
    setup_initial_particles(&particles, &weights, N_features, xtrue);
    setup_initial_Q_R();  // modifies global variables

    int *ftag;
    int *da_table;
    setup_landmarks(&ftag, &da_table, N_features);

    Vector2d *z;  // This is a dynamic array of Vector2d - see https://stackoverflow.com/a/13597383/5616591
    Vector2d *zf;
    Vector2d *zn;
    int *idf, *ftag_visible;
    setup_measurements(&z, &zf, &zn, &idf, &ftag_visible, N_features);

//    if ( SWITCH_PREDICT_NOISE ) {
//        printf("Sampling from predict noise usually OFF for FastSLAM 2.0\n");	
//    }
 
    if ( SWITCH_SEED_RANDOM ) {
        srand( SWITCH_SEED_RANDOM );
    }	

    double dt        = DT_CONTROLS; // change in time btw predicts
    double dtsum     = 0;           // change in time since last observation
    double T         = 0;
    int iwp          = 0;           // index to first waypoint
    double G         = 0;           // initialize steering angle
    double V         = V_; 
    size_t Nf_visible = 0;

    // Main loop
    while ( iwp != -1 ) {

        //////////////////////////////////////////////////////////////////
        // Prediction
        //////////////////////////////////////////////////////////////////

        predict_update(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);

        /////////////////////////////////////////////////////////////////


        //Update time
        dtsum = dtsum + dt;
        T+=dt;

        // Observation condition
        if ( dtsum >= DT_OBSERVE ) {
            dtsum = 0;
            
            //////////////////////////////////////////////////////////////
            // Observation
            //////////////////////////////////////////////////////////////

            observe_update(lm, N_features, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);

            //////////////////////////////////////////////////////////////
        }
    }

    cleanup_landmarks(&ftag, &da_table);
    cleanup_measurements(&z, &zf, &zn, &idf, &ftag_visible);
    *particles_ = particles;
    *weights_ = weights;
}
