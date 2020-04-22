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

void fastslam1_sim( double* lm, const size_t lm_rows, const size_t lm_cols, 
                    double* wp, const size_t wp_rows, const size_t wp_cols, 
                    Particle *particle ) 
{
    if (SWITCH_PREDICT_NOISE) {
        printf("Sampling from predict noise usually OFF for FastSLAM 2.0\n");	
    }

    // normally initialized configfile.h
    // Matrix2d declared in configfile.h
    Q[0] = pow(sigmaV,2);
    Q[1] = 0;
    Q[2] = 0;
    Q[3] = pow(sigmaG,2); 

    // Matrix2d declared in configfile.h
    R[0] = sigmaR*sigmaR;
    R[1] = 0;
    R[2] = 0; 
    R[3] = sigmaB*sigmaB;

    double veh[2][3] = { {0,-WHEELBASE,-WHEELBASE}, {0,-1,1} };

    // vector of particles (their count will change)
    Particle *particles = (Particle *) malloc( NPARTICLES * sizeof(Particle) );
    for (size_t i = 0; i < NPARTICLES; i++) {
        initParticle(&particle[i], lm_cols); // Nf = lm_cols; is this true ?
    }

    // initialize particle weights as uniform
    double uniformw = 1.0 / (double) NPARTICLES; 
    double* w = (double *) malloc( NPARTICLES * sizeof(double) );
    for (size_t i = 0; i < NPARTICLES; i++) {
        w[i] = uniformw;
        particles[i].w = &w[i];
    }

    double xtrue[3] = {0.0, 0.0, 0.0};

    double dt = DT_CONTROLS; // change in time btw predicts
    double dtsum = 0; // change in time since last observation

    size_t *ftag = (size_t*) malloc( lm_cols * sizeof(size_t) ); // identifier for each lm
    for (size_t i = 0; i < lm_cols; i++) {
        ftag[i] = i; 
    }

    // data association table
    double *da_table = (double *) malloc( lm_cols * sizeof(double) );
    for (size_t i = 0; i < lm_cols; i++) {
        da_table[i] = -1;
    }

    int iwp = 0; // index to first waypoint
    double G = 0; // initial steer angle
    // MatrixXd plines; // will later change to list of points

    if (SWITCH_SEED_RANDOM !=0) {
        srand(SWITCH_SEED_RANDOM);
    } 		

    Matrix2d Qe;
    Matrix2d Re;
    copy(Q, 4, Qe);
    copy(R, 4, Re);

    if (SWITCH_INFLATE_NOISE ==1) {
        scal(Qe, 4, 2.0, Qe);
        scal(Re, 4, 2.0, Re);
    }

    int* ftag_visible;
    double *z; // vector<Vector2d> -> range and bearings of visible landmarks

    // Main loop
    while ( iwp != -1 ) {
        compute_steering(xtrue, wp, N_wp, AT_WAYPOINT, RATEG, MAXG, dt, &iwp, &G);
        if (iwp == -1 && NUMBER_LOOPS > 1) {
            iwp = 0;
            NUMBER_LOOPS = NUMBER_LOOPS - 1;
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
            /* if (SWITCH_HEADING_KNOWN) { */
            /* for (int j=0; j< particles[i].xf().size(); j++) { */
            /* Vector2d xf_j = particles[i].xf()[j]; */
            /* xf_j[2] = xtrue[2]; */
            /* particles[i].setXfi(j,xf_j); */	
            /* } */           
            /* } */
        }

        //Observe step
        dtsum = dtsum+dt;
        if (dtsum >= DT_OBSERVE) {
            dtsum=0;

            //Compute true data, then add noise
            ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	

            //z is the range and bearing of the observed landmark
            z = get_observations(xtrue,lm,ftag_visible,MAX_RANGE);
            add_observation_noise(z,R,SWITCH_SENSOR_NOISE);

            if (!z.empty()){
                plines = make_laser_lines(z,xtrue);
            }

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
        }
    }
    cout<<"done with all functions and will return particles"<<endl<<flush;
    return particles;
}

//rb is measurements
//xv is robot pose
/*MatrixXd make_laser_lines(vector<Vector2d> rb, Vector3d xv) 
{
    if (rb.empty()) {
        return MatrixXd(0,0);
    }

    int len = rb.size();
    MatrixXd lnes(4,len);

    MatrixXd globalMat(2,rb.size());
    int j;
    for (j=0; j<globalMat.cols();j++) {
        globalMat(0,j) = rb[j][0]*cos(rb[j][1]); 
        globalMat(1,j) = rb[j][0]*sin(rb[j][1]);   	
    }

    TransformToGlobal(globalMat,xv);

    for (int c=0; c<lnes.cols();c++) {
        lnes(0,c) = xv(0);
        lnes(1,c) = xv(1);
        lnes(2,c) = globalMat(0,c);
        lnes(3,c) = globalMat(1,c);
    }

    return line_plot_conversion(lnes);
}*/
