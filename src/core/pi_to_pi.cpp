#include "pi_to_pi.h"

#include <math.h>

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base Implementation, unit test, Benchmark added
 *       While and if implementation -> probably doesn't get much faster
 * ToDo: Maybe a unrolled wihle implementation? Debatable though...
 ****************************************************************************/

double pi_to_pi(double ang) {
    return pi_to_pi_while(ang);
}

/*****************************************************************************
 * PERFORMANCE STATUS BASE
 * Work, best: 1 negate + 2 mult + 4 fl-comp = 7 flops
 * Work, worst: 2 negate + 6 mult + 1 div + 1 floor + 4 fl-comp + 2 add=16 fl
 * Memory moved: 0 (register)
 * Cycles: 32 cyc
 * Performance: 0.18 flops / cyc
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

double pi_to_pi_base_flops(double ang) {
    double flop_count = 2*tp.mul + 2*tp.negation * 4*tp.doublecomp;
    if ((ang <= (-2* M_PI)) || (ang > (2*M_PI))) {
        int n=floor(ang/(2*M_PI));
        ang = ang-n*(2*M_PI);    
        flop_count+= tp.floor + tp.div + tp.mul;
    }
    if (ang > M_PI) {
        ang = ang - (2*M_PI);
        flop_count+= tp.add + tp .mul;
    }
    if (ang <= -M_PI) {
        ang = ang + (2*M_PI);
        flop_count+= tp.add + tp .mul;
    }
    return flop_count;
}

double pi_to_pi_base_memory(double ang) {
    return 0.0;
}

double pi_to_pi_active_flops(double ang) {
    return pi_to_pi_base_flops(ang);
}

double pi_to_pi_active_memory(double ang) {
    return 0.0;
}

void pi_to_pi_arr(double* angles,const size_t n) 
{
    for (int i=0; i<n; i++) {
        angles[i] = pi_to_pi(angles[i]);
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

/*****************************************************************************
 * PERFORMANCE STATUS BASE
 * Work, best: 2 fl-comp = 2 flops
 * Work, worst: 4 fl-comp + 1 fl-sub = 5 flops if angle in [-2PI, 2PI], otherwise potentially arbitrarily long
 * Memory moved: 0 (register)
 * Cycles: 20 cyc
 * Performance: 0.29 flops / cyc
 * Optimal: I think to
 * Status: Pretty good
 ****************************************************************************/
double pi_to_pi_while(double ang) {
    if (ang > M_PI) {
        while (ang > M_PI) ang -= 2*M_PI;
    }
    else if (ang <= -M_PI) {
        while (ang <= -M_PI) ang += 2*M_PI;
    }
    return ang;
}