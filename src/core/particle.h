#pragma once
#include "typedefs.h"
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>


typedef struct Particle {
	double* w; //! importance weight
	double xv[3]; //! robot pose: x,y,theta (heading dir)
	double Pv[9]; //! control inputs, i.e. velocities
	int Nf; //! Max Number of features (you should hardly ever need this)
	int Nfa; //! Actual number of known features (you should hardly ever need this)
	double* xf; //! 2d means of EKF
	double* Pf; //! covariance matrices for EKF 
	void (*del) (Particle& p);
	void (*delPtr) (Particle* p);
	void (*set_xv) (Particle* p, Vector3d xv);
	void (*set_Pv) (Particle* p, Matrix3d Pv);
	void (*set_xfi) (Particle* p, Vector3d xf, int index);
	void (*set_Pfi) (Particle* p, Matrix3d Pf, int index);
	int index;
} Particle;

/*!
    Initialize existing particle.
    @param[out] p   Pointer to a particle.
	@param[in]  Nf  	Maximum number of expected features.
 */
void initParticle(Particle* p, const size_t Nf);

/*! Deep copy contents of one particle into the other */
void copyParticle(const Particle& p_ref, Particle& p_target);

/*!
    Initialize existing particle with particle_index.
    @param[out] p   Pointer to a particle.
	@param[in]  Nf  	Maximum number of expected features.
	@param[in]  index  	Index for the constructed particle.
 */
void initParticle(Particle* p, const size_t Nf, int particle_index);

/*!
    Constructor for an empty particle.
    @param[in]  Nf  Maximum number of expected features.
	@return 	p	Constructed particle.
 */
Particle* newParticle(const size_t Nf);

/*!
    Constructor for an empty particle with index.
    @param[in]  Nf  	Maximum number of expected features.
	@param[in]  index  	Index for the constructed particle.
	@return 	p	Constructed particle.
 */
Particle* newParticle(const size_t Nf, int particle_index);

/*!
    Free particle.
    @param[in]  p   Ref to a particle.
 */
void delParticle (Particle& p);

/*!
    Free particle behind ptr.
    @param[in]  p   Pointer to a particle.
 */
void delParticlePtr (Particle* p);

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

