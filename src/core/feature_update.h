#pragma once

#include "KF_cholesky_update.h"
#include "compute_jacobians.h"
#include "particle.h"
#include "pi_to_pi.h"

/*! Updates the state of a particle given a list of measurements.
 *
 *  [Runtime depends on pi_to_pi and KF_cholesky, otherwise mostly memory
 *  operations. Switch to mask!]
 *      @param[out] 	particle 	Particle to be updated.
 *      @param[in] 		z		    list of measurements conditioned on the particle.
 *      @param[out] 	idf 	    Index of known landmarks.
 *      @param[out] 	R	 	    Covariance Matrix of measurements.
 */
void feature_update(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]);

void feature_update_base(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]);


void feature_update_active(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]);

// Work / Memory instrumenting
double feature_update_base_flops(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]);

double feature_update_base_memory(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]);

double feature_update_active_flops(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]);

double feature_update_active_memory(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]);
