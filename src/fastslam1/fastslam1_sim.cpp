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
#include "fastrand.h"

void fastslam1_sim( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_)
{
    fastslam1_sim_active(lm, lm_rows, lm_cols, wp, wp_rows, wp_cols, particles_, weights_);
}

double fastslam1_sim_base_flops( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_) {
    const size_t N_features = lm_rows;
    const size_t N_waypoints = wp_rows;
    NUMBER_LOOPS=_NUMBER_LOOPS;
    Particle *particles;
    double *weights;

    double weights_copy[NPARTICLES];
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

    double flop_count = 0;

//    if ( SWITCH_PREDICT_NOISE ) {
//        printf("Sampling from predict noise usually OFF for FastSLAM 2.0\n");	
//    }
 
    srand( SWITCH_SEED_RANDOM );
#ifdef __AVX2__
    avx_xorshift128plus_init(1,1);
#endif

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
        flop_count+= predict_update_base_flops(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);
        predict_update_base(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);

        /////////////////////////////////////////////////////////////////


        //Update time
        dtsum = dtsum + dt;
        T+=dt;

        flop_count+= 2* tp.add;

        // Observation condition
        if ( dtsum >= DT_OBSERVE ) {
            dtsum = 0;
            
            //////////////////////////////////////////////////////////////
            // Observation
            //////////////////////////////////////////////////////////////
            for (size_t i = 0; i < N_features; i++) {
                ftag_visible[i] = ftag[i];
            }

            flop_count += get_observations_base_flops(xtrue, MAX_RANGE, lm, N_features, ftag_visible, &Nf_visible, z);
            get_observations_base(xtrue, MAX_RANGE, lm, N_features, ftag_visible, &Nf_visible, z); // Nf_visible = number of visible features
            
            
            //z is the range and bearing of the observed landmark

            flop_count+= observe_update_base_flops(lm, N_features, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);
            observe_update_base(lm, N_features, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);

            //////////////////////////////////////////////////////////////
        }
    }

    cleanup_landmarks(&ftag, &da_table);
    cleanup_measurements(&z, &zf, &zn, &idf, &ftag_visible);
    cleanup_particles(&particles, &weights);
    return flop_count;
}

double fastslam1_sim_base_memory( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_){
    const size_t N_features = lm_rows;
    const size_t N_waypoints = wp_rows;
    NUMBER_LOOPS=_NUMBER_LOOPS;
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

    double memory_moved = 0.0;

//    if ( SWITCH_PREDICT_NOISE ) {
//        printf("Sampling from predict noise usually OFF for FastSLAM 2.0\n");	
//    }
 
    srand( SWITCH_SEED_RANDOM );
#ifdef __AVX2__
    avx_xorshift128plus_init(1,1);
#endif

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
        memory_moved+= predict_update_base_memory(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);
        predict_update_base(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);

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
            for (size_t i = 0; i < N_features; i++) {
                ftag_visible[i] = ftag[i];
            }

            //z is the range and bearing of the observed landmark
            memory_moved+= get_observations_base_memory(xtrue, MAX_RANGE, lm, N_features, ftag_visible, &Nf_visible, z);
            get_observations_base(xtrue, MAX_RANGE, lm, N_features, ftag_visible, &Nf_visible, z); // Nf_visible = number of visible features
            

            memory_moved+= observe_update_base_memory(lm, N_features, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);
            observe_update_base(lm, N_features, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);

            //////////////////////////////////////////////////////////////
        }
    }

    cleanup_landmarks(&ftag, &da_table);
    cleanup_measurements(&z, &zf, &zn, &idf, &ftag_visible);
    cleanup_particles(&particles, &weights);
    return memory_moved;
}


double fastslam1_sim_active_flops( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_){
    const size_t N_features = lm_rows;
    const size_t N_waypoints = wp_rows;
    NUMBER_LOOPS=_NUMBER_LOOPS;
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

    double flop_count = 0.0;

//    if ( SWITCH_PREDICT_NOISE ) {
//        printf("Sampling from predict noise usually OFF for FastSLAM 2.0\n");	
//    }
 
    srand( SWITCH_SEED_RANDOM );
#ifdef __AVX2__
    avx_xorshift128plus_init(1,1);
#endif

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
        flop_count+= predict_update_active_flops(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);
        predict_update_base(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);

        /////////////////////////////////////////////////////////////////


        //Update time
        dtsum = dtsum + dt;
        T+=dt;

        flop_count+= 2* tp.add;

        // Observation condition
        if ( dtsum >= DT_OBSERVE ) {
            dtsum = 0;
            
            //////////////////////////////////////////////////////////////
            // Observation
            //////////////////////////////////////////////////////////////
            for (size_t i = 0; i < N_features; i++) {
                ftag_visible[i] = ftag[i];
            }

            //z is the range and bearing of the observed landmark
            flop_count+= get_observations_base_flops(xtrue, MAX_RANGE, lm, N_features, ftag_visible, &Nf_visible, z);
            get_observations_base(xtrue, MAX_RANGE, lm, N_features, ftag_visible, &Nf_visible, z); // Nf_visible = number of visible features

            flop_count+= observe_update_flops(lm, N_features, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);
            observe_update_base(lm, N_features, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);
            	
            //////////////////////////////////////////////////////////////
        }
    }

    cleanup_landmarks(&ftag, &da_table);
    cleanup_measurements(&z, &zf, &zn, &idf, &ftag_visible);
    cleanup_particles(&particles, &weights);
    return flop_count;
}

double fastslam1_sim_active_memory( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_) {
    const size_t N_features = lm_rows;
    const size_t N_waypoints = wp_rows;
    NUMBER_LOOPS = _NUMBER_LOOPS;
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

    double memory_moved = 0.0;

//    if ( SWITCH_PREDICT_NOISE ) {
//        printf("Sampling from predict noise usually OFF for FastSLAM 2.0\n");	
//    }
 
    srand( SWITCH_SEED_RANDOM );
#ifdef __AVX2__
    avx_xorshift128plus_init(1,1);
#endif

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
        memory_moved+= predict_update_active_memory(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);
        predict_update_base(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);

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
            for (size_t i = 0; i < N_features; i++) {
                ftag_visible[i] = ftag[i];
            }

            //z is the range and bearing of the observed landmark
            memory_moved += get_observations_base_memory(xtrue, MAX_RANGE, lm, N_features, ftag_visible, &Nf_visible, z);
            get_observations_base(xtrue, MAX_RANGE, lm, N_features, ftag_visible, &Nf_visible, z); // Nf_visible = number of visible features

            memory_moved += observe_update_base_memory(lm, N_features, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);
            observe_update_base(lm, N_features, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);

            //////////////////////////////////////////////////////////////
        }
    }

    cleanup_landmarks(&ftag, &da_table);
    cleanup_measurements(&z, &zf, &zn, &idf, &ftag_visible);
    cleanup_particles(&particles, &weights);
    return memory_moved;
}


void fastslam1_sim_base( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_)
{
    const size_t N_features = lm_rows;
    const size_t N_waypoints = wp_rows;
    NUMBER_LOOPS=_NUMBER_LOOPS;
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
 
    srand( SWITCH_SEED_RANDOM );
#ifdef __AVX2__
    avx_xorshift128plus_init(1,1);
#endif

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

        predict_update_base(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);

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
            for (size_t i = 0; i < N_features; i++) {
                ftag_visible[i] = ftag[i];
            }

            //z is the range and bearing of the observed landmark
    
            get_observations_base(xtrue, MAX_RANGE, lm, N_features, ftag_visible, &Nf_visible, z); // Nf_visible = number of visible features

            observe_update_base(lm, N_features, xtrue, *R, ftag, 
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

void fastslam1_sim_active( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle **particles_, double** weights_)
{
    const size_t N_features = lm_rows;
    const size_t N_waypoints = wp_rows;
    NUMBER_LOOPS=_NUMBER_LOOPS;
    Particle *particles;
    double *weights;
    // double *xv; We get them from the configfile
    // double *Pv; We get them from the configfile

    Vector3d xtrue   = {0,0,0};
    setup_initial_particles_and_pose(&particles, &weights, &xv, &Pv, NPARTICLES, N_features, xtrue);
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
    srand( SWITCH_SEED_RANDOM );
#ifdef __AVX2__
    avx_xorshift128plus_init(1,1);
#endif

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

            for (size_t i = 0; i < N_features; i++) {
                ftag_visible[i] = ftag[i];
            }

            //z is the range and bearing of the observed landmark
    
            get_observations_base(xtrue, MAX_RANGE, lm, N_features, ftag_visible, &Nf_visible, z); // Nf_visible = number of visible features

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

void fastslam1_sim_base_VP(double* lm, const size_t lm_rows, const size_t lm_cols, 
        const size_t N_features, Particle **particles_, double** weights_)
{
    const size_t N_waypoints = 0;
    double *wp = NULL; // dummy

    Particle *particles;
    double *weights;
    Vector3d xtrue = {-67.6493, -41.7142, 35.5*M_PI/180};
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

    srand( SWITCH_SEED_RANDOM );
#ifdef __AVX2__
    avx_xorshift128plus_init(1,1);
#endif

    double dt        = DT_CONTROLS; // change in time btw predicts
    double dtsum     = 0;           // change in time since last observation
    double T         = lm[0];
    int iwp          = 0;           // index to first waypoint
    double G         = 0;           // initialize steering angle
    double V         = 0.0;//V_; 
    size_t Nf_visible = 0;

    // Main loop
    int index = 0;
    while ( index < lm_rows -1 ) {

        //////////////////////////////////////////////////////////////////
        // Prediction
        //////////////////////////////////////////////////////////////////

        dt = (lm[4*(index+1)] - lm[4*(index)]) / 1000.0;

        V = lm[4*(index) +2];
        G = lm[4*(index) +3];
        predict_update_VP_base(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);

        /////////////////////////////////////////////////////////////////

        //Update time
        dtsum = dtsum + dt;
        T += dt*1000.0;

        // Observation condition
        if ( lm[4*(index+1) + 1] > -1 ) {
            dtsum = 0;

            ///Setup z, ftag_visible, Nf_visible
            Nf_visible = 0;
            while ( lm[4*(index+1) + 1] > -1 ) {
                index++;
                double r = lm[4*(index) +2];
                double phi = lm[4*(index) +3];

                z[Nf_visible][0] = r;
                z[Nf_visible][1] = pi_to_pi_base(phi - M_PI_2);

                ftag_visible[Nf_visible] = lm[4*(index) +1] -1;
                Nf_visible++;
            }

            //////////////////////////////////////////////////////////////
            // Observation
            //////////////////////////////////////////////////////////////

            observe_update_base(lm, N_features, xtrue, *R, ftag, 
                    da_table, ftag_visible, z, &Nf_visible, zf, idf, 
                    zn, particles, weights);

            dt = (lm[4*(index+1)] - lm[4*(index)]) / 1000.0;

            T += dt*1000.0;

            predict_update_VP_base(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G, particles);

            //////////////////////////////////////////////////////////////
        }
        index++;
    }

    cleanup_landmarks(&ftag, &da_table);
    cleanup_measurements(&z, &zf, &zn, &idf, &ftag_visible);
    *particles_ = particles;
    *weights_ = weights;
}

void fastslam1_sim_active_VP(double* lm, const size_t lm_rows, const size_t lm_cols, 
        const size_t N_features, Particle **particles_, double** weights_)
{
    const size_t N_waypoints = 0;
    double *wp = NULL; // dummy

    Particle *particles;
    double *weights;
    Vector3d xtrue = {-67.6493, -41.7142, 35.5*M_PI/180};
    setup_initial_particles_and_pose(&particles, &weights, &xv, &Pv, NPARTICLES, N_features, xtrue);
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

    srand( SWITCH_SEED_RANDOM );
#ifdef __AVX2__
    avx_xorshift128plus_init(1,1);
#endif

    double dt        = DT_CONTROLS; // change in time btw predicts
    double dtsum     = 0;           // change in time since last observation
    double T         = lm[0];
    int iwp          = 0;           // index to first waypoint
    double G         = 0;           // initialize steering angle
    double V         = 0.0;//V_; 
    size_t Nf_visible = 0;

    // Main loop
    int index = 0;
    while ( index < lm_rows -1 ) {

        //////////////////////////////////////////////////////////////////
        // Prediction
        //////////////////////////////////////////////////////////////////

        dt = (lm[4*(index+1)] - lm[4*(index)]) / 1000.0;

        V = lm[4*(index) +2];
        G = lm[4*(index) +3];
        predict_update_VP_active(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);

        /////////////////////////////////////////////////////////////////

        //Update time
        dtsum = dtsum + dt;
        T += dt*1000.0;

        // Observation condition
        if ( lm[4*(index+1) + 1] > -1 ) {
            dtsum = 0;

            ///Setup z, ftag_visible, Nf_visible
            Nf_visible = 0;
            while ( lm[4*(index+1) + 1] > -1 ) {
                index++;
                double r = lm[4*(index) +2];
                double phi = lm[4*(index) +3];

                z[Nf_visible][0] = r;
                z[Nf_visible][1] = pi_to_pi_base(phi - M_PI_2);

                ftag_visible[Nf_visible] = lm[4*(index) +1] -1;
                Nf_visible++;
            }

            //////////////////////////////////////////////////////////////
            // Observation
            //////////////////////////////////////////////////////////////

            observe_update(lm, N_features, xtrue, *R, ftag, 
                    da_table, ftag_visible, z, &Nf_visible, zf, idf, 
                    zn, particles, weights);

            dt = (lm[4*(index+1)] - lm[4*(index)]) / 1000.0;

            T += dt*1000.0;

            predict_update_VP_active(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G, particles);

            //////////////////////////////////////////////////////////////
        }
        index++;
    }

    cleanup_landmarks(&ftag, &da_table);
    cleanup_measurements(&z, &zf, &zn, &idf, &ftag_visible);
    *particles_ = particles;
    *weights_ = weights;
}
