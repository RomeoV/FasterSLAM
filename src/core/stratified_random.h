#pragma once
#include "typedefs.h"

/*!
    Stratified random vector (Sampling with subpopulations).
    @param[in]  N     Interval size (inverse).
    @param[out] di    Generated array of random stratified numbers.
 */
void stratified_random(int N, double* di);

/*!
    Generates a random, uniformly sampled number between [0,1].
    @return random uniform number in [0,1].
 */
double unifRand();

