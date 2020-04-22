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
// #include "line_plot_conversion.h" //don't need this?

int counter = 0;

typedef struct {
    double xtrue[3];
    double veh[2][3];
    double V;
    double alpha; // initial steer angle, used to be called G in yglee
} Vehicle;

void setup_initial_particles(Particle **p_, double **w_, const size_t N_features);
void setup_initial_Q_R();
void setup_initial_vehicle(Vehicle* v);
void setup_landmarks(int **ftag_, double **da_table_, const size_t N_features);
void setup_measurements(double **z, int **ftag_visible, const size_t N_features);
void destroy_particles(Particle **p_, double **w_);
void destroy_landmarks(int **ftag_, double **da_table_);
void destroy_measurements(double **z, int **ftag_visible);

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

    if ( SWITCH_PREDICT_NOISE ) {
        printf("Sampling from predict noise usually OFF for FastSLAM 2.0\n");	
    }
 
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
    std::cout << "done with all functions and has updated particles"<<std::endl<<std::flush;
    // return particles;
}

void setup_initial_particles(Particle **p_, double **w_, const size_t N_features) {
    *p_ = (Particle *) malloc( NPARTICLES * sizeof(Particle) );
    Particle *p = *p_;
    for (size_t i = 0; i < NPARTICLES; i++) {
        initParticle(&p[i], N_features);
    }
    const double uniform_w = 1.0 / NPARTICLES; 
    *w_ = (double *) malloc( NPARTICLES * sizeof(double) );
    double *w = *w_;
    for (size_t i = 0; i < NPARTICLES; i++) {
        w[i] = uniform_w;
        p[i].w = &w[i];
    }
}

void destroy_particles(Particle** p_, double** w_) {
    Particle* p = *p_;
    for (size_t i = 0; i < NPARTICLES; i++) {
        delParticleMembers(&p[i]);
    }
    free(p_);
    free(w_);
}

void setup_initial_Q_R() {
    // Matrix2d Q, R are declared in configfile.h
    Q[0][0] = pow(sigmaV, 2);
    Q[0][1] = 0;
    Q[1][0] = 0;
    Q[1][1] = pow(sigmaG, 2); 
    
    R[0][0] = sigmaR * sigmaR;
    R[0][1] = 0;
    R[1][0] = 0;
    R[1][1] = sigmaB * sigmaB;

    if ( SWITCH_INFLATE_NOISE == 1 ) {
        scal(*Q, 4, 2.0, *Q);
        scal(*R, 4, 2.0, *R);
    }
}

void setup_initial_vehicle(Vehicle* v) {
    Vector3d xtrue_initial = {0.0, 0.0, 0.0};
    copy(v->xtrue, 3, xtrue_initial);

    double veh_initial[2][3] = { {0,-WHEELBASE,-WHEELBASE}, {0,-1,1} };
    copy(*(v->veh), 2*3, *veh_initial);

    v->V = 0; // \todo
    v->alpha = 0.0;  // along global coordinate x axis
}

void setup_landmarks(int **ftag_, double **da_table_, const size_t N_features) {
    *ftag_ = (int*) malloc( N_features * sizeof(int) ); // identifier for each lm
    *da_table_ = (double *) malloc( N_features * sizeof(double) ); // data association table
    for (size_t i = 0; i < N_features; i++) {
        (*ftag_)[i] = i; 
        (*da_table_)[i] = -1;
    }
}

void destroy_landmarks(int **ftag_, double **da_table_) {
    free(*ftag_);
    free(*da_table_);
}

void setup_measurements(double **z, int **ftag_visible, const size_t N_features) {
    *z = (double*) malloc( 2 * N_features * sizeof(double) );
    *ftag_visible = (int*) malloc( N_features * sizeof(int) );
}

void destroy_measurements(double **z, size_t **ftag_visible) {
    free(*z);
    free(*ftag_visible);
}
