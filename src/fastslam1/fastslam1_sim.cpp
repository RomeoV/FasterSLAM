#include <iostream>
#include <math.h>
#include "fastslam1_sim.h"
#include "add_control_noise.h" // don't have it yet
#include "get_observations.h"
#include "add_observation_noise.h"
#include "TransformToGlobal.h" // don't have it yet
// #include "line_plot_conversion.h" //don't have it yet
#include "data_associate_known.h" //don't have it yet
#include "feature_update.h"
#include "resample_particles.h"
#include "add_feature.h"
#include "compute_weight.h"
#include "predict.h"

int counter = 0;

typedef struct {
    double xtrue[3];
    double veh[2][3];
    double alpha; // initial steer angle, used to be called G in yglee
} Vehicle;

// allocs and inits array of particles (its size will change)
// allocs and inits particle weights as uniform
void setup_initial_particles(Particle **p_, double **w_);
void setup_initial_Q_R();
void setup_initial_vehicle(Vehicle* v);
void setup_landmarks(size_t **ftag_, double **da_table_, const size_t N_features);
void setup_measurements(double **z, size_t **ftag_visible);
void destroy_particles(Particle **p_, double **w_);
void destroy_landmarks(size_t **ftag_, double **da_table_);
void destroy_measurements(double **z, size_t **ftag_visible);

void fastslam1_sim( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle *particle ) 
{
    const size_t N_features = lm_rows;

    Particle *particles;
    double *weights;
    setup_initial_particles(&particles, &weights);

    setup_inital_Q_R();  // modifies global variables

    Vehicle vehicle_gt;
    setup_initial_vehicle(&vehicle_gt);

    size_t *ftag;
    double *da_table;
    setup_landmarks(&ftag, &da_table);

    double *z; // range and bearings of visible landmarks ( vector<Vector2d> )
    size_t *ftag_visible;
    setup_measurements(&z, &ftag_visible);

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
        
        compute_steering(xtrue, wp, N_wp, AT_WAYPOINT, RATEG, MAXG, dt, &iwp, &G);
        if ( iwp == -1 && NUMBER_LOOPS > 1 ) {
            iwp = 0;
            NUMBER_LOOPS--;
        }
        predict_true(V, G, WHEELBASE, dt, xv);

        // add process noise
        double* VnGn = new double[2];        
        add_control_noise(V, G, Q, SWITCH_CONTROL_NOISE, VnGn); // TODO
        double Vn = VnGn[0];
        double Gn = VnGn[1];

        // Predict step	
        for (unsigned int i = 0; i < NPARTICLES; i++) {
            predict(particles[i], Vn, Gn, Qe, WHEELBASE, dt, SWITCH_PREDICT_NOISE);
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
            ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	

            //z is the range and bearing of the observed landmark
            z = get_observations(xtrue,lm,ftag_visible,MAX_RANGE);
            add_observation_noise(z,R,SWITCH_SENSOR_NOISE);

//            if (!z.empty()){
//                plines = make_laser_lines(z,xtrue);
//            }

            //Compute (known) data associations
            int Nf = particles[0].xf().size();
            vector<int> idf;
            vector<Vector2d> zf;
            vector<Vector2d> zn;            

            bool testflag= false;
            data_associate_known(z,ftag_visible,da_table,Nf,zf,idf,zn);

            // perform update
            for (int i =0; i<NPARTICLES; i++) {
                if (!zf.empty()) { //observe map features
                    double w = compute_weight(particles[i],zf,idf,R);
                    w = particles[i].w()*w;
                    particles[i].setW(w);
                    feature_update(particles[i],zf,idf,R);
                }
                if (!zn.empty()) {
                    add_feature(particles[i], zn, R);
                }
            }

            resample_particles(particles,NEFFECTIVE,SWITCH_RESAMPLE);

            if (VnGn) { 
                delete[] VnGn;
            }
///////////////////////////////////////////////////////////////////////////////////
        }
    }
    cout << "done with all functions and will return particles"<<endl<<flush;
    return particles;
}

void setup_initial_particles(Particle **p_, double **w_) {
    *p_ = (Particle *) malloc( NPARTICLES * sizeof(Particle) );
    Particle *p = *p_;
    for (size_t i = 0; i < NPARTICLES; i++) {
        initParticle(&p[i], lm_rows);
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
        delParticleMembers(p[i]);
    }
    free(p_);
    free(w_);
}

void setup_initial_Q_R() {

    // Matrix2d Q, R are declared in configfile.h
    Q[0] = pow(sigmaV, 2);
    Q[1] = 0;
    Q[2] = 0;
    Q[3] = pow(sigmaG, 2); 
    
    R[0] = sigmaR * sigmaR;
    R[1] = 0;
    R[2] = 0;
    R[3] = sigmaB * sigmaB;

    if ( SWITCH_INFLATE_NOISE == 1 ) {
        scal(Q, 4, 2.0, Q);
        scal(R, 4, 2.0, R);
    }
}

void setup_initial_vehicle(Vehicle* v) {
    Vector3d xtrue_initial = {0.0, 0.0, 0.0};
    copy(v.xtrue, 3, xtrue_initial);

    double veh_initial[2][3] = { {0,-WHEELBASE,-WHEELBASE}, {0,-1,1} };
    copy(v.veh, 2*3, veh_initial);

    v.alpha = 0.0;  // along global coordinate x axis
}

void setup_landmarks(size_t **ftag_, double **da_table_, const size_t N_features) {
    *ftag = (size_t*) malloc( N_features * sizeof(size_t) ); // identifier for each lm
    *da_table = (double *) malloc( N_features * sizeof(double) ); // data association table
    for (size_t i = 0; i < N_features; i++) {
        (*ftag)[i] = i; 
        (*da_table)[i] = -1;
    }
}

void destroy_landmarks(size_t **ftag_, double **da_table_) {
    free(*ftag_);
    free(*da_table);
}

void setup_measurements(double **z, size_t **ftag_visible) {
    *z = (double*) malloc( 2 * N_features * sizeof(double) );
    *ftag_visible = (size_t*) malloc( N_features * sizeof(double) );
}

void destroy_measurements(double **z, size_t **ftag_visible) {
    free(*z);
    free(*ftag_visible);
}