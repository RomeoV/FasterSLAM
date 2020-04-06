#ifndef ADD_OBSERVATION_NOISE_H
#define ADD_OBSERVATION_NOISE_H

#include <Eigen/Dense>
#include <iostream>
#include <vector>

using namespace Eigen;
using namespace std;

namespace nRandMat{
	MatrixXd randn(int m, int n); //Gaussian distribution
	MatrixXd rand(int m, int n); //Standard random
}

/*!
    Adds random observation noise to the observation vector. [Compute intensive]
    @param[out] z        	Landmark measurements / observations [meter, radians].
    @param[in]  R        	Covariance matrix of observation (diagonal).
    @param[in]  addnoise	Flag if obersvation noise should be added.
 */
void add_observation_noise(vector<Vector2d> &z, Matrix2d R, int addnoise);

#endif //ADD_OBSERVATION_NOISE_H
