#include "observe_update.h"

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

#include "configfile.h"

#include <stddef.h>

void observe_update(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    observe_update_base(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf,
                        zn, particles, weights);
}

void observe_update_base(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

    //z is the range and bearing of the observed landmark
    
    get_observations_base(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    
    add_observation_noise_base(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible -> idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known_base(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions

    // perform update
    for (size_t i = 0; i < NPARTICLES; i++) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            double w = compute_weight_base(&particles[i], zf, count_zf, idf, R);
            w *= weights[i];
            weights[i] = w;
            feature_update_base(&particles[i], zf, idf, count_zf, R);
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature_base(&particles[i], zn, count_zn, R);
        }
    }

    resample_particles_base(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}

void observe_update_active(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    // Compute true data, then add noise
    // ftag_visible = vector<int>(ftag); //modify the copy, not the ftag	
    // memcpy(ftag_visible, ftag, N_features*sizeof(int));
    for (size_t i = 0; i < N_features; i++) {
        ftag_visible[i] = ftag[i];
    }

    //z is the range and bearing of the observed landmark
    
    get_observations(xtrue, MAX_RANGE, lm, N_features, ftag_visible, Nf_visible, z); // Nf_visible = number of visible features
    //print(*z,*Nf_visible,2,std::cout);
    if ( *Nf_visible == 0 ) {
        return;
    }
    
    add_observation_noise(z, *Nf_visible, R, SWITCH_SENSOR_NOISE);

    //Compute (known) data associations
    const int Nf_known = particles[0].Nfa; // >= Nf_visible -> idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions

    // perform update
    for (size_t i = 0; i < NPARTICLES; i++) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            double w = compute_weight(&particles[i], zf, count_zf, idf, R);
            w *= weights[i];
            weights[i] = w;
            feature_update(&particles[i], zf, idf, count_zf, R);
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(&particles[i], zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}