#ifndef KF_CHOLESKY_UPDATE_H
#define KF_CHOLESKY_UPDATE_H

#include <Eigen/Dense>

using namespace Eigen;


/*!
    Update Kalman Filter using Cholesky decomposition. [Some LinAlg, Nik could do that ;)]
	@param[out] x 	    State vector.
    @param[out] P 	    Covariance matrix of state.
	@param[in]  v   	velocity at current state.
    @param[in]  R    	Covariance matrix of measurements.
    @param[in]  H    	Jacobian of h wrt feature states.
 */
void KF_cholesky_update(Vector2d &x,Matrix2d &P,Vector2d v,Matrix2d R,Matrix2d H);

#endif
