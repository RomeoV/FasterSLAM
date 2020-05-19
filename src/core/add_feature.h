#pragma once
#include "particle.h"
#include "typedefs.h"

/*! Adds a feature observation to a particle.

    [Memory-heavy, leave for later, but switch to mask]
    @param[out] particle Particle where the feature is added to.
    @param[in]  z        Landmark measurements / observations in robot
   coordinates [meter, radians].
    @param[in]  R        Covariance matrix of observation noises -> configfile.h
 */
void add_feature(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R);

void add_feature_base(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R);

void add_feature_active(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R);

// Work / Memory instrumenting
double add_feature_base_flops(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R);

double add_feature_base_memory(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R);

double add_feature_active_flops(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R);

double add_feature_active_memory(Particle* particle, Vector2d z[], size_t N_z, Matrix2d R);
