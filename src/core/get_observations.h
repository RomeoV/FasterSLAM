#pragma once
#include "typedefs.h"

/******************************************************************
 * All Compute-heavy
 ************************************************************/
/*!
    Updates the state of a particle given a list of measurements.
	@param[in] 	x 	    State vector (x,y,angle).
	@param[in] 	lm		Matrix of landmark positions. 2 x NumLM
	@param[out] idf 	Index of known landmarks. (updated)
	@param[in] 	rmax 	Maximal distance between x and lm(:,i) to be considered visible.

    @return Vector of distance and bearing to each landmark.
 */
vector<Vector2d> get_observations(Vector3d x, MatrixXd lm, vector<int> &idf, double rmax);

/*!
    Computes which landmarks are visible by comparing the distance x<->lm(:,i) to rmax.
    Only returns visible landmarks in lm and idf.
	@param[in] 	x 	    State vector (x,y,angle).
	@param[out] lm		Matrix of landmark positions. 2 x NumLM
	@param[out] idf 	Index of known landmarks. (updated)
	@param[in] 	rmax 	Maximal distance between x and lm(:,i) to be considered visible.
 */
void get_visible_landmarks(Vector3d x, MatrixXd &lm, vector<int> &idf, double rmax);

/*!
    Computes distance and bearing between x and all landmarks in lm.
	@param[in] 	x 	    State vector (x,y,angle).
	@param[out] lm		Matrix of landmark positions. 2 x NumLM

    @return Vector  of 2d vectors (distance, angle) between x and each landmark in lm. 
 */
vector<Vector2d> compute_range_bearing(Vector3d x, MatrixXd lm);

/*!
    Finds all visible landmarks given current state.
    @param[in] 	dx      Vector of x-distances to all landmarks.
    @param[in] 	dy      Vector of x-distances to all landmarks.
    @param[in] 	size    Number of landmarks.
    @param[in]  phi     Angle of current state.
    @param[in] 	rmax    Maximal distance between x and lm(:,i) to be considered visible.

    @return Vector landmark indices which are visible and its size.
 */
//! index should be preallocated with size equal to size
void find2(const double *dx, const double *dy, const size_t size, 
           const double phi, const double rmax, double *index, size_t *index_size);

