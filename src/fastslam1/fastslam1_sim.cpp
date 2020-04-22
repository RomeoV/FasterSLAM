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
#include "transform_to_global.h"
#include "fastslam1_utils.h"
// #include "line_plot_conversion.h" //don't need this?

int counter = 0;

void fastslam1_sim( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle *particle ) 
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

    double *z; // range and bearings of visible landmarks ( vector<Vector2d> )
    int *ftag_visible;
    setup_measurements(&z, &ftag_visible, N_features);

//    if ( SWITCH_PREDICT_NOISE ) {
//        printf("Sampling from predict noise usually OFF for FastSLAM 2.0\n");	
//    }
 
    if ( SWITCH_SEED_RANDOM ) {
        srand( SWITCH_SEED_RANDOM );
    }	

    double dt        = DT_CONTROLS; // change in time btw predicts
    double dtsum     = 0;           // change in time since last observation
    int iwp          = 0;           // index to first waypoint

    // Main loop
    while ( iwp != -1 ) {
//////////////////////////////////////////////////////////////////////////
// ground_truth_update()
//////////////////////////////////////////////////////////////////////////
        
        compute_steering(vehicle_gt.xtrue, wp, N_waypoints, AT_WAYPOINT, RATEG, MAXG, dt, &iwp, &vehicle_gt.alpha);
        if ( iwp == -1 && NUMBER_LOOPS > 1 ) {
            iwp = 0;
            NUMBER_LOOPS--;
        }
        Vector2d xv;
        predict_true(vehicle_gt.V, vehicle_gt.alpha, WHEELBASE, dt, xv);

        // add process noise
        double* VnGn = new double[2];        
        add_control_noise(vehicle_gt.V, vehicle_gt.alpha, *Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
        double Vn = VnGn[0];
        double Gn = VnGn[1];

        // Predict step	
        for (unsigned int i = 0; i < NPARTICLES; i++) {
            predict(&particles[i], Vn, Gn, *Q, dt);
        }

//////////////////////////////////////////////////////////////////////////
        //Observe step
        dtsum = dtsum+dt;
        if (dtsum >= DT_OBSERVE) {
            dtsum=0;
///////////////////////////////////////////////////////////////////////////////////
// observe()
//////////////////////////////////////////////////////////////////////////
            //Compute true data, then add noise
            // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
            memcpy(ftag_visible, ftag, N_features*sizeof(int));

            //z is the range and bearing of the observed landmark
            size_t N_measurements;
            get_observations(vehicle_gt.xtrue, MAX_RANGE, lm, N_features, &ftag_visible, &N_measurements, z);
            add_observation_noise(z, N_measurements, *R,SWITCH_SENSOR_NOISE);

//            if (!z.empty()){
//                plines = make_laser_lines(z,xtrue);
//            }

            //Compute (known) data associations
            int Nf = particles[0].Nfa;
            vector<int> idf;
            vector<Vector2d> zf;
            vector<Vector2d> zn;            

            bool testflag = false;
            data_associate_known(z, ftag_visible, da_table, Nf, zf, idf, zn);

            // perform update
            for (int i = 0; i < NPARTICLES; i++) {
                if ( !zf.empty() ) { //observe map features
                    double w = compute_weight(&particles[i], zf, idf, *R);
                    w = (*particles[i].w)*w;
                    *particles[i].w = w;
                    feature_update(&particles[i], zf, idf, *R);
                }
                if ( !zn.empty() ) {
                    add_feature(&particles[i], zn, *R);
                }
            }

            // resample_particles(particles,NEFFECTIVE); \todo what's with NEFFECTIVE?
            resample_particles(particles, NPARTICLES, weights);

            if (VnGn) { 
                delete[] VnGn;
            }
///////////////////////////////////////////////////////////////////////////////////
        }
    }
}
