




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

#include "read_input_file.h"

#include "estimate_json.cpp"
#include "ground_truth_json.cpp"

#include "json.hpp"

int main (int argc, char *argv[])
{
	std::string output_filename = "robot_trace.json";
    std::string ground_truth_filename = "ground_truth.json";

    double *lm; // landmark positions
	double *wp; // way points
    size_t lm_rows, wp_rows;

	if (argc < 2)
		return -1;

	read_input_file(argv[1], &lm, &wp, lm_rows, wp_rows);

    const size_t N_features = lm_rows;
    const size_t N_waypoints = wp_rows;

    Particle *particles;
    double *weights;
    Vector3d xtrue   = {0,0,0};

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
 
    if ( SWITCH_SEED_RANDOM ) {
        srand( SWITCH_SEED_RANDOM );
    }	
    uint64_t init_state[8] = {1,1,1,1,1,1,1,1};
    avx_xorshift128plus_init(1,1);
    uint64_t init_seq[8] = {1,3,5,7,9,11,13,15};
    pcg32_srand(1,1);
    avx2_pcg32_srand(init_state, init_seq);

    double dt        = DT_CONTROLS; // change in time btw predicts
    double dtsum     = 0;           // change in time since last observation
    double T         = 0;
    int iwp          = 0;           // index to first waypoint
    double G         = 0;           // initialize steering angle
    double V         = V_; 
    size_t Nf_visible = 0;

    // Simulation parameters
    bool observe=true;
	int id=0;
    int observe_id =0;
    nlohmann::json particle_trace = {{"timesteps", nlohmann::json::array()}};
	nlohmann::json ground_truth = {{"timesteps", nlohmann::json::array()}};
	ground_truth.update(ground_truth_keypoints_json(wp, lm, N_waypoints, N_features));

    // Main loop
    while ( iwp != -1 ) {

        //Optional
        if (id == 60000){
            break;
        }


        //JSON
        if (observe) {
            particle_trace["timesteps"] += estimate_step_json(particles, weights, NPARTICLES, T, 
            ftag_visible, Nf_visible, da_table, N_features, id) ;
            observe=false;
            observe_id++;
        }

        ground_truth["timesteps"]+= ground_truth_step_json(xtrue, T, id, iwp, G);
        id++;
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
            observe=true;
            
            //////////////////////////////////////////////////////////////
            // Observation
            //////////////////////////////////////////////////////////////

            observe_update(lm, N_features, xtrue, *R, ftag, 
            da_table, ftag_visible, z, &Nf_visible, zf, idf, 
            zn, particles, weights);

            //////////////////////////////////////////////////////////////
        }
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
	free(wp);
}
