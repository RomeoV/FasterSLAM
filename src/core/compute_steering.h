#pragma once
#include "typedefs.h"
#include <stddef.h>
/*!
    Motion Model along given set of waypoints.
    Required for prediction step. [Simple computations, could start now]
    @param[in] 	x 		State vector. Pointer to array of length 3. (x,y,angle) [EDIT]
    @param[in] 	wp 		waypoints. Pointer to 1D-array of length 2*n. [EDIT]
    @param[in]  N_wp    Number of waypoints. [NEW]
    @param[in] 	minD 	minimum distance to current waypoint before switching to next.
    @param[in] 	rateG 	max steering rate (rad/s).
    @param[in] 	maxG 	max steering angle (rad).
    @param[in] 	dt 		timestep.
    @param[out] iwp 	index to current waypoint. ToDo: Handle iwp=-1 (second round)
    @param[out] G 		current steering angle.
 */
void compute_steering(cVector3d x, double* wp, const size_t N_wp, const double minD, 
                      const double rateG, const double maxG, const double dt, 
                      int* iwp, double* G);

void compute_steering_base(cVector3d x, double* wp, const size_t N_wp, const double minD, 
                      const double rateG, const double maxG, const double dt, 
                      int* iwp, double* G);

                      
// Utils
FlopCount compute_steering_base_flops(cVector3d x, double* wp, const size_t N_wp,
                                   const double minD, const double rateG,
                                   const double maxG, const double dt,
                                   int* iwp, double* G); 

FlopCount compute_steering_active_flops(cVector3d x, double* wp, const size_t N_wp,
                                   const double minD, const double rateG,
                                   const double maxG, const double dt,
                                   int* iwp, double* G); 

double compute_steering_base_memory(cVector3d x, double* wp, const size_t N_wp,
                                    const double minD, const double rateG,
                                    const double maxG, const double dt,
                                    int* iwp, double* G); 

double compute_steering_active_memory(cVector3d x, double* wp, const size_t N_wp,
                                    const double minD, const double rateG,
                                    const double maxG, const double dt,
                                    int* iwp, double* G); 
