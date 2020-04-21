#include <iostream>
#include <math.h>
#include <vector>

#include "fastslam1_sim.h"
#include "core/add_control_noise.h"
#include "core/get_observations.h"
#include "core/add_observation_noise.h"
#include "core/TransformToGlobal.h"
#include "core/line_plot_conversion.h"
#include "core/data_associate_known.h"
#include "core/feature_update.h"
#include "core/resample_particles.h"
#include "core/add_feature.h"
#include "compute_weight.h"
#include "predict.h"

using namespace config;
using namespace std;

int counter=0;
vector<Particle> fastslam1_sim(MatrixXd lm, MatrixXd wp) 
{
    if (SWITCH_PREDICT_NOISE) {
	printf("Sampling from predict noise usually OFF for FastSLAM 2.0\n");	
    }

    //normally initialized configfile.h
    Q << pow(sigmaV,2), 0,
      0 , pow(sigmaG,2);

    R << sigmaR*sigmaR, 0, 
      0, sigmaB*sigmaB;

    double veh[2][3] = {{0,-WHEELBASE,-WHEELBASE},{0,-1,1}};

    //vector of particles (their count will change)
    vector<Particle> particles(NPARTICLES);
    for (int i=0; i<particles.size(); i++){
	particles[i] = Particle();
    }

    //initialize particle weights as uniform
    double uniformw = 1.0/(double)NPARTICLES;    
    for (unsigned int p = 0; p < NPARTICLES; p++) {
	particles[p].setW(uniformw);
    }

    Vector3d xtrue(3);
    xtrue.setZero();

    double dt = DT_CONTROLS; //change in time btw predicts
    double dtsum = 0; //change in time since last observation

    vector<int> ftag; //identifier for each landmark
    for (int i=0; i< lm.cols(); i++) {
	ftag.push_back(i); 
    }

    //data ssociation table
    VectorXd da_table(lm.cols());
    for (int s=0; s<da_table.size(); s++) {
	da_table[s] = -1;
    }

    int iwp = 0; //index to first waypoint
    double G = 0; //initial steer angle
    MatrixXd plines; //will later change to list of points

    if (SWITCH_SEED_RANDOM !=0) {
	srand(SWITCH_SEED_RANDOM);
    } 		

    Matrix2d Qe = Matrix2d(Q);
    Matrix2d Re = Matrix2d(R);

    if (SWITCH_INFLATE_NOISE ==1) {
	Qe = 2*Q;
	Re = 2*R;
    }

    vector<int> ftag_visible;
    vector<Vector2d> z; //range and bearings of visible landmarks

    //Main loop
    while (iwp !=-1) {
	compute_steering(xtrue, wp, iwp, AT_WAYPOINT, G, RATEG, MAXG, dt);
	if (iwp ==-1 && NUMBER_LOOPS > 1) {
	    iwp = 0;
	    NUMBER_LOOPS = NUMBER_LOOPS-1;
	}
	predict_true(xtrue,V,G,WHEELBASE,dt);

	//add process noise
	double* VnGn = new double[2];        
	add_control_noise(V,G,Q,SWITCH_CONTROL_NOISE,VnGn);
	double Vn = VnGn[0];
	double Gn = VnGn[1];

	//Predict step	
	for (unsigned int i=0; i< NPARTICLES; i++) {
	    predict(particles[i],Vn,Gn,Qe,WHEELBASE,dt,SWITCH_PREDICT_NOISE);
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
	    
	    //perform update
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
MatrixXd make_laser_lines(vector<Vector2d> rb, Vector3d xv) 
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
}
