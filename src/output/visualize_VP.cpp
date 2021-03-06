




#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <iomanip>

#include "linalg.h"
#include "fastrand.h"
#include "predict_update.h"
#include "observe_update.h"
#include "configfile.h"
#include "fastslam1_utils.h"
#include "pi_to_pi.h"

#include "read_input_file.h"

#include "estimate_json.cpp"
#include "ground_truth_json.cpp"

#include "json.hpp"

int main (int argc, char *argv[])
{
	std::string output_filename = "robot_trace_old.json";
    std::string ground_truth_filename = "ground_truth_old.json";

    double *lm; // landmark positions
	double *wp; // way points
    size_t lm_rows, wp_rows;

	if (argc < 2)
		return -1;
    NPARTICLES = 1000;
	//read_input_file(argv[1], &lm, &wp, lm_rows, wp_rows);

	
	size_t N_features = 0;
    read_sequential_input_file(argv[1], &lm, lm_rows, N_features);

	const size_t lm_cols = 4;
    const size_t N_waypoints = 0;

    Particle *particles;
    double *weights;
    Vector3d xtrue   =  {-67.6493, -41.7142, 35.5*M_PI/180};

    if (argc > 2) {
        setup_initial_particles(&particles, &weights, N_features, xtrue);
    } else {
        setup_initial_particles_and_pose(&particles, &weights, &xv, &Pv, NPARTICLES, N_features, xtrue);
    }

    setup_initial_Q_R();  // modifies global variables

    int *ftag;
    int *da_table;
    setup_landmarks(&ftag, &da_table, N_features);

    Vector2d *z; 
    Vector2d *zf;
    Vector2d *zn;
    int *idf, *ftag_visible;
    setup_measurements(&z, &zf, &zn, &idf, &ftag_visible, N_features);

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

    // Simulation parameters
    bool observe=true;
	int id=0;
    int observe_id =0;

	int index = 0;

    nlohmann::json particle_trace = {{"timesteps", nlohmann::json::array()}};
	nlohmann::json ground_truth = {{"timesteps", nlohmann::json::array()}};
	ground_truth.update(ground_truth_keypoints_json(wp, lm, N_waypoints, 0));


    // Main loop
    while ( index < lm_rows -1 ) {
        // Optional
        // if (index > 20000) {
        //     break;
        // }


        //JSON
        if (observe) {
            particle_trace["timesteps"] += estimate_step_json(particles, weights, NPARTICLES, T, 
            ftag_visible, Nf_visible, da_table, N_features, id) ;
            observe=false;
            observe_id++;
        }

        ground_truth["timesteps"]+= ground_truth_step_json(xtrue, V, id, index, G);
        id++;
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
        T+=dt*1000.0;

        // Observation condition
        if ( lm[4*(index+1) + 1] > -1) {
            dtsum = 0;
            observe=true;


			///Setup z, ftag_visible, Nf_visible
			Nf_visible = 0;
			while (lm[4*(index+1) + 1] > -1) {
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

			T+=dt*1000.0;

			predict_update_VP_base(wp, N_waypoints, V, *Q, dt, NPARTICLES, xtrue, &iwp, &G,particles);

            //////////////////////////////////////////////////////////////
        }
		index++;
    }

    // Write JSON
    std::ofstream of(output_filename);
    of << particle_trace;

	std::ofstream of_gt(ground_truth_filename);
    of_gt << ground_truth;

    cleanup_landmarks(&ftag, &da_table);
    cleanup_measurements(&z, &zf, &zn, &idf, &ftag_visible);

    if (argc > 2) {
        cleanup_particles(&particles, &weights);
    } else {
        cleanup_particles_and_pose(&particles, &weights, &xv, &Pv, NPARTICLES);
    }
    

	free(lm);
	//free(wp);
}
