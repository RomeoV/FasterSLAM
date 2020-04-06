#pragma once

/*!
    Adds random observation noise to the observation vector. [Compute intensive]
    @param[out] z        	Landmark measurements / observations [meter, radians].
    @param[in]  R        	Covariance matrix of observation (diagonal).
    @param[in]  addnoise	Flag if obersvation noise should be added.
 */
void add_observation_noise(double *z, const int lenz, cMatrix2d R, const int addnoise);
