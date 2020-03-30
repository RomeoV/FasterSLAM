#pragma once
#include "typedefs.h"

/*!
    Stratified resampling based on weights w. Interval size is N_w. [If we assume w is
    a fixed array, where the array is passed as an argument, we can start working here.
    keep should be a mask then.]
    @param[in]   w          Array of weights.
    @param[in]   N_w        Size of w and keep array.
    @param[out]  Neff       = 1/sum(w_i^2).
    @param[out]  keep       Array of ints that indicate whether w[i] should be kept or not.
 */
void stratified_resample(double* w, int N_w, double* Neff, int* keep);

/*!
    Returns an array where all weights leq its index are summed up. 
    w_out(i) = sum_j=0->j=i-1(w_in(j))
    @param[out] w       Array of weights.
    @param[in]  N_w   Length of array.
 */
void cumsum(double* w, int N_w);
