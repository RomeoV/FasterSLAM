#pragma once
#include "typedefs.h"

/*!
    Clips ang  to range [-pi,pi].
    @param[out] angle   angle in radiants.
 */
double pi_to_pi(double ang);
double pi_to_pi_base(double ang);


/*!
    Clips all angles in angle (radiants) to range [-pi,pi] (in place!)
    @param[out] angles   Vector of angles (radiants).
    @param[in]  n       Length of angles array. [NEW]
 */
void pi_to_pi_arr(double* angles, const size_t n);
void pi_to_pi_arr_base(double* angles, const size_t n);

/*! This is supposed to be a faster version of pi_to_pi using fmod */
double pi_to_pi_fmod(double ang);

/*! This is a version that only works on (-2PI, 2PI]. It is even a bit faster than pi_to_pi_while(double) though*/
double pi_to_pi_nongeneral(double ang);

/*! This is a version that checks in a while loop if the angle is smaller or greater and keeps +/-ing 2PI */
double pi_to_pi_while(double ang);

double pi_to_pi_base_flops(double ang);
double pi_to_pi_base_memory(double ang);

double pi_to_pi_active_flops(double ang);
double pi_to_pi_active_memory(double ang);
