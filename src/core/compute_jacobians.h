#pragma once
#include "typedefs.h"
#include "particle.h"


/*! Thrun03 Eq. 38, 39 and 60
 *
 *  Note that the sensor model is defined as follows:
 *  $ g(\theta, s_t) = \begin{bmatrix} \sqrt{(x_{\theta} - x_v)^2 + (y_{\theta} - y_v)^2}  \\
 *                                     \arctan(\frac{y_{\theta} - y_v}{x_{\theta} - x_v}) \end{bmatrix} $
 * 
 *  First computes all predicted observations in relative coordinates to the vehicle.
 *  Then computes the jacobians given a particle state and predict observations. [Compute-Intensive, Switch to mask]
 *  @param[in]   Particle   Particle for which the jacobian should be computed.
 *  @param[in]   idf        Feature indices.
 *  @param[in]   R          Covariance matrix of observation (diagonal).
 *  @param[out]  zp         vector of predicted observation (given the new vehicle state)
 *  @param[out]  Hv         Jacobian of h wrt vehicle states
 *  @param[out]  Hf         Jacobian of h wrt feature states
 *  @param[out]  Sf         Measurement covariance of feature observation given the vehicle.
 */
void compute_jacobians(Particle* particle, 
        int idf[], 
        size_t N_z,
        Matrix2d R, 
        Vector2d* zp, //measurement (range, bearing)
        Matrix23d* Hv, // jacobians of function h (deriv of h wrt pose)
        Matrix2d* Hf, // jacobians of function h (deriv of h wrt mean)
        Matrix2d* Sf) //measurement covariance
;
