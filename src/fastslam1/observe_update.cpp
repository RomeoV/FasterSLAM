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
#include <math.h>

#include "configfile.h"

#include <stddef.h>

void observe_update(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) {
    observe_update_active(lm, N_features, xtrue, R, ftag, da_table, ftag_visible, z, Nf_visible, zf, idf,
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
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known_base(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions

    // perform update
    for (size_t i = 0; i < NPARTICLES; i++) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            // std::cout<<"HiBase"<<std::endl;
            // TODO: precompute the jacobian of compute_weight_base and feature_update_base
            // this will half the processing need of compute_jacobians
            Vector2d zp[count_zf] __attribute__ ((aligned(32)));
            Matrix23d Hv[count_zf] __attribute__ ((aligned(32)));
            Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
            Matrix2d Sf[count_zf] __attribute__ ((aligned(32)));
            double w = compute_weight_base(&particles[i], zf, count_zf, idf, R, zp, Hv, Hf, Sf);
            w *= weights[i];
            weights[i] = w;
            feature_update_base(&particles[i], zf, idf, count_zf, R, zp, Hv, Hf, Sf);
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
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions

    // perform update
    for (size_t i = 0; i < NPARTICLES; i++) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            Vector2d zp[count_zf] __attribute__ ((aligned(32)));
            Matrix23d Hv[count_zf] __attribute__ ((aligned(32)));
            Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
            Matrix2d Sf[count_zf] __attribute__ ((aligned(32)));
            compute_jacobians_fast(&particles[i], idf, count_zf, R, zp, Hv, Hf, Sf);
            double w = compute_weight_active(&particles[i], zf, count_zf, idf, R, zp, Hv, Hf, Sf);
            w *= weights[i];
            weights[i] = w;
            feature_update_active(&particles[i], zf, idf, count_zf, R, zp, Hv, Hf, Sf);
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(&particles[i], zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}

void observe_update_fast(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
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
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions
    Vector2d zp[count_zf] __attribute__ ((aligned(32)));
    Matrix23d Hv[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf[count_zf] __attribute__ ((aligned(32)));

    Vector2d v[count_zf]  __attribute__ ((aligned(32)));
    Vector2d feat_diff[count_zf];  // difference btw feature prediciton and
                                        // measurement (used to update mean)
    // perform update

    //double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;
    double den, num;
    Matrix2d S_inv;
    Vector2d S_inv_v;
    double vT_S_inv_v;
    for (size_t i = 0; i < NPARTICLES; i++) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            // std::cout<<"HiFast"<<std::endl;
            //COMPUTE JACOBIANS
            compute_jacobians_fast(particles + i, idf, count_zf, R, zp, Hv, Hf, Sf);
            //END COMPUTE_JACOBIANS
            double w = 1.0;
            for (size_t j = 0; j < count_zf; j++) {
                sub(z[j], zp[j], 2, feat_diff[j]);
                feat_diff[j][1] = pi_to_pi(feat_diff[j][1]);
                inv_2x2(Sf[j], S_inv);
                
                mv_2x2(S_inv, feat_diff[j], S_inv_v);
                mul(feat_diff[j], S_inv_v, 1, 2, 1, &vT_S_inv_v); // TODO in linalg

                den = 2 * M_PI * sqrt(determinant_2x2(Sf[j]));
                num = exp(-0.5 * vT_S_inv_v);
                w *= (double)num / (double)den;

                
                KF_cholesky_update(particles[i].xf + 2 * idf[j], particles[i].Pf + 4 * idf[j], 
                                feat_diff[j], R, Hf[j]);
                set_xfi(particles+i, particles[i].xf + 2 * idf[j], idf[j]);
                set_Pfi(particles+i, particles[i].Pf + 4 * idf[j], idf[j]);
            }
            weights[i]*=w;
            
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(particles+i, zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}

void observe_update_simplified(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
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
    const int Nf_known = particles[0].Nfa; // >= Nf_visible . idz_size
    size_t count_zf = 0;
    size_t count_zn = 0;
    data_associate_known(z, ftag_visible, *Nf_visible, da_table, Nf_known, zf, idf, &count_zf, zn, &count_zn); // TODO Rewrite/fix bugs + create test for this functions
    Vector2d zp[count_zf] __attribute__ ((aligned(32)));
    Matrix23d Hv[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Hf[count_zf] __attribute__ ((aligned(32)));
    Matrix2d Sf[count_zf] __attribute__ ((aligned(32)));

    Vector2d v[count_zf]  __attribute__ ((aligned(32)));
    Vector2d feat_diff[count_zf];  // difference btw feature prediciton and
                                        // measurement (used to update mean)
    // perform update

    //double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;
    double den, num;
    Matrix2d S, ST, S_inv;
    Vector2d S_inv_v;
    double vT_S_inv_v;
    for (size_t i = 0; i < NPARTICLES; i++) {
        if ( count_zf != 0 ) { //observe map features ( !zf.empty() )
            // std::cout<<"HiFast"<<std::endl;
            //COMPUTE JACOBIANS
            compute_jacobians_fast(particles + i, idf, count_zf, R, zp, Hv, Hf, Sf);
            //END COMPUTE_JACOBIANS
            double w = 1.0;
            for (size_t j = 0; j < count_zf; j++) {
                sub(z[j], zp[j], 2, feat_diff[j]);
                feat_diff[j][1] = pi_to_pi(feat_diff[j][1]);
                inv_2x2(Sf[j], S_inv);
                
                mv_2x2(S_inv, feat_diff[j], S_inv_v);
                mul(feat_diff[j], S_inv_v, 1, 2, 1, &vT_S_inv_v); // TODO in linalg

                den = 2 * M_PI * sqrt(determinant_2x2(Sf[j]));
                num = exp(-0.5 * vT_S_inv_v);
                w *= (double)num / (double)den;

                
                KF_cholesky_update(particles[i].xf + 2 * idf[j], particles[i].Pf + 4 * idf[j], 
                                feat_diff[j], R, 
                                Hf[j]);
                set_xfi(particles+i, particles[i].xf + 2 * idf[j], idf[j]);
                set_Pfi(particles+i, particles[i].Pf + 4 * idf[j], idf[j]);
            }
            weights[i]*=w;
            
        }
        if ( count_zn != 0 ) { // !zn.empty() 
            add_feature(particles+i, zn, count_zn, R);
        }
    }

    resample_particles(particles, NPARTICLES, weights, NEFFECTIVE, SWITCH_RESAMPLE);            
}


