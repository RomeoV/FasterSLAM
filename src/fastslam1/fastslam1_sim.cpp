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

int counter = 0;

void fastslam1_sim( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particle ) 
{
    const size_t N_features = lm_rows;
    const size_t N_waypoints = wp_rows;

    Particle *particles;
    double *weights;
    setup_initial_particles(&particles, &weights, N_features);
    setup_initial_Q_R();  // modifies global variables
    Vehicle vehicle_gt;
    setup_initial_vehicle(&vehicle_gt);

    int *ftag;
    double *da_table;
    setup_landmarks(&ftag, &da_table, N_features);

    Vector2d (*z);  // This is a dynamic array of Vector2d - see https://stackoverflow.com/a/13597383/5616591
    Vector2d (*zf);
    Vector2d (*zn);
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
    int iwp          = 0;           // index to first waypoint
    double G         = 0;           // initialize steering angle

    // Main loop
    while ( iwp != -1 ) {
//////////////////////////////////////////////////////////////////////////
// ground_truth_update()
//////////////////////////////////////////////////////////////////////////

        compute_steering(vehicle_gt.xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, &iwp, &G);
        if ( iwp == -1 && NUMBER_LOOPS > 1 ) {
            iwp = 0;
            NUMBER_LOOPS--;
        }
        predict_true(vehicle_gt.V, vehicle_gt.xtrue[2], WHEELBASE, dt, vehicle_gt.xtrue);

        // add process noise
        double VnGn[2];
        add_control_noise(vehicle_gt.V, vehicle_gt.xtrue[2], *Q, SWITCH_CONTROL_NOISE, VnGn); // TODO

        // Predict step	
        for (size_t i = 0; i < NPARTICLES; i++) {
            predict(&particles[i], VnGn[0], VnGn[1], *Q, dt);
        }

//////////////////////////////////////////////////////////////////////////
        //Observe step
        dtsum = dtsum + dt;
        if ( dtsum >= DT_OBSERVE ) {
            dtsum = 0;
///////////////////////////////////////////////////////////////////////////////////
// observe()
//////////////////////////////////////////////////////////////////////////
            // Compute true data, then add noise
            // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
            // memcpy(ftag_visible, ftag, N_features*sizeof(int));
            for (size_t i = 0; i < N_features; i++) {
                ftag_visible[i] = ftag[i];
            }

            //z is the range and bearing of the observed landmark
            size_t N_measurements;
            get_observations(vehicle_gt.xtrue, MAX_RANGE, lm, N_features, ftag_visible, &N_measurements, z); // N_measurements = number of visible features
            
            if ( N_measurements == 0 ) {
                continue;
            }
            
            add_observation_noise(z, N_measurements, *R, SWITCH_SENSOR_NOISE);

            //Compute (known) data associations
            const int Nf_known = particles[0].Nfa; // >= N_measurements -> idz_size
            size_t count_zf = 0;
            size_t count_zn = 0;
            data_associate_known(z, ftag_visible, N_measurements, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions

            // perform update
            for (size_t i = 0; i < NPARTICLES; i++) {
                if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
                    double w = compute_weight(&particles[i], zf, count_zf, idf, *R);
                    w *= *( particles[i].w );
                    *( particles[i].w ) = w;
                    feature_update(&particles[i], zf, idf, count_zf, *R);
                }
                if ( count_zn != 0 ) { // !zn.empty() 
                    add_feature(&particles[i], zn, count_zn, *R);
                }
            }

            resample_particles(particles, NPARTICLES, weights);

///////////////////////////////////////////////////////////////////////////////////
        }
    }
    *particle = particles;
}
