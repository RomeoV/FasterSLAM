#include "pi_to_pi.h"

#include <math.h>


/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base Implementation, unit test, Benchmark added
 * ToDo: Start optimizing
 ****************************************************************************/

double pi_to_pi(double ang) {
    return pi_to_pi_base(ang);
}

/*****************************************************************************
 * PERFORMANCE STATUS BASE
 * Work, best: 1 negate + 2 mult + 4 fl-comp = 7 flops
 * Work, worst: 2 negate + 6 mult + 1 div + 1 floor + 4 fl-comp + 2 add=16 fl
 * Memory moved: 0 (register)
 * Cycles: 42 cyc
 * Performance: 0.33 flops / cyc
 * Optimal: TBD
 * Status: Baseline
 ****************************************************************************/

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

double pi_to_pi_nongeneral(double ang) {
    if (ang > M_PI) ang -= 2*M_PI;
    else if (ang <= M_PI) ang += 2*M_PI;
    return ang;
}

double pi_to_pi_while(double ang) {
    if (ang > M_PI) {
        while (ang > M_PI) ang -= 2*M_PI;
    }
    else if (ang <= M_PI) {
        while (ang <= M_PI) ang += 2*M_PI;
    }
    return ang;
}