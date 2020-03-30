#pragma once
#include "typedefs.h"

/*!
    Clips ang  to range [-pi,pi].
    @param[out] angle   angle in radiants.
 */
double pi_to_pi(double ang);


/*!
    Clips all angles in angle (radiants) to range [-pi,pi].
    @param[out] angles   Vector of angles (radiants).
    @param[in]  n       Length of angles array. [NEW]
 */
void pi_to_pi_arr(double* angles, int n);


