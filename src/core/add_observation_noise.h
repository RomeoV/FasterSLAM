#pragma once
#include "typedefs.h"

/*!
  Adds random observation noise to the observation vector. [Compute intensive]
  @param[out] z        	Landmark measurements / observations [meter, radians].
  @param[in]  R        	Covariance matrix of observation (diagonal).
  @param[in]  addnoise	Flag if obersvation noise should be added.
  */
void add_observation_noise(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise);
double add_observation_noise_flops(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise);
double add_observation_noise_memory(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise);

void add_observation_noise_base(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise);
double add_observation_noise_base_flops(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise);
double add_observation_noise_base_memory(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise);
