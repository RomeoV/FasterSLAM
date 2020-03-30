#pragma once
#include "typedefs.h"

/*!
    Update Kalman Filter using Cholesky decomposition. 
    @param[out] x 	    State vector.
    @param[out] P 	    Covariance matrix of state.
	@param[in]  v   	Velocity at current state.
    @param[in]  R    	Covariance matrix of measurements.
    @param[in]  H    	Jacobian of h wrt feature states.
 */
void KF_cholesky_update(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);
