#pragma once
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "typedefs.h"

typedef struct Particle {
  double* w;    //! importance weight
  double* xv;  //! robot pose: x,y,theta (heading dir)
  double* Pv;  //! control inputs, i.e. velocities
  int Nf;    //! Max Number of features (you should hardly ever need this)
  int Nfa;   //! Actual number of known features (you should hardly ever need
                //! this)
  double* xf;   //! 2d means of EKF in cartesian world coordinates
  double* Pf;   //! covariance matrices for EKF (=Extended Kalman Filter) in polar robot coordinates
} Particle;


/*! Deep copy contents of one particle into the other */
void copyParticle(const Particle* p_ref, Particle* p_target);

/*!
    Initialize existing particle with particle_index. xv & Pv point to external memory.
    @param[out] p   Pointer to a particle.
        @param[in]  Nf  	        Maximum number of expected features.
        @param[in]  xv_initial      Array[3] of initial state vector.
 */
void initParticle_prealloc(Particle* p, const int Nf, const double* xv_initial);


/*!
    Initialize existing particle with particle_index. Alloc xv and Pv.
    @param[out] p   Pointer to a particle.
        @param[in]  Nf  	        Maximum number of expected features.
        @param[in]  xv_initial      Array[3] of initial state vector.
 */
void initParticle(Particle* p, const int Nf, const double* xv_initial);

/*!
    Constructor for an empty particle with index. Allocs xv and Pv.
    @param[out] p   Pointer to a particle.
        @param[in]  Nf  	        Maximum number of expected features.
        @param[in]  xv_initial      Array[3] of initial state vector.
 */
Particle* newParticle(const int Nf, const double* xv_initial);

/*!
    Free particle.
    @param[in]  p   Pointer to a particle.
 */
void delParticleMembers(Particle* p);

/*!
    Free particle.
    @param[in]  p   Pointer to a particle.
 */
void delParticleMembers_prealloc(Particle* p);

/*!
    Free particle behind ptr.
    @param[in]  p   Pointer to a particle.
 */
void delParticleMembersAndFreePtr(Particle* p);

/*!
    Free particle behind ptr.
    @param[in]  p   Pointer to a particle.
 */
void delParticleMembersAndFreePtr_prealloc (Particle* p);

/*!
    Set covariance matrix of this particles state vector.
    @param[out] p  	Pointer to a particle.
        @param[in]	Pv	Covariane matrix of state vector.
 */
void set_Pv(Particle* p, Matrix3d Pv);

/*!
    Set state vector of this particle.
    @param[out] p  	Pointer to a particle.
        @param[in]	xv	State vector.
 */
void set_xv(Particle* p, Vector3d xv);

/*!
    Set feature vector of this particle.
    @param[out] p  		Pointer to a particle.
        @param[in]	xf		State vector of respective feature.
        @param[in]	index	Index where the feature should be set.
 */
void set_xfi(Particle* p, Vector2d xf, int index);

/*!
    Set feature covariance matrix of this particle and feature at index.
    @param[out] p  		Pointer to a particle.
        @param[in]	Pf		Covariane matrix of respective feature.
        @param[in]	index	Index where the feature should be set.
 */
void set_Pfi(Particle* p, Matrix2d Pf, int index);
