#include "pi_to_pi.h"

#include <math.h>


/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base Implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: TBD
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/

double pi_to_pi(double ang) {
    return pi_to_pi_base(ang);
}

double pi_to_pi_base(double ang) 
{
    if ((ang <= (-2* M_PI)) || (ang > (2*M_PI))) {
        int n=floor(ang/(2*M_PI));
        ang = ang-n*(2*M_PI);    
    }
    if (ang > M_PI) {
        ang = ang - (2*M_PI);
    }
    if (ang <= -M_PI) {
        ang = ang + (2*M_PI);
    }
    return ang;
}

void pi_to_pi_arr(double* angles,const size_t n) 
{
    pi_to_pi_arr_base(angles, n);
}

void pi_to_pi_arr_base(double* angles,const size_t n) 
{
    for (int i=0; i<n; i++) {
        angles[i] = pi_to_pi_base(angles[i]);
    }
}

double pi_to_pi_fmod(double ang) {
    /* I think this can be done more efficiently with
     * ```
     * theta = ((theta+pi)%(2*pi)-pi)
     * ```
     * Note that `a%b` is the remainder, not the modulo, so we have to add `b` if we are less than zero.
     * More info here:
     * https://stackoverflow.com/questions/13683563/whats-the-difference-between-mod-and-remainder
     */
    double tmp = fmod(ang+M_PI, 2*M_PI);
    if (tmp <= 0) tmp += 2*M_PI;
    return tmp - M_PI;
}
