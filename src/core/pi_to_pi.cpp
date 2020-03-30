#include "pi_to_pi.h"

#include <math.h>
#include <stdio.h>


/*****************************************************************************
 * OPTIMIZATION STATUS
 * Last Worked on: 30.03.2020
 * Done: Base Implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: TBD
 * Memory moved: TBD
 * Cycles: Not measured.
 * Performance: Not measured.
 * Optimal: Not measured.
 * Status: Not started.
 ****************************************************************************/

double pi_to_pi(double ang) 
{
    int n;
    if ((ang < (-M_PI)) || (ang > (M_PI))) {
        n=floor(ang/(2*M_PI));
        ang = ang-n*(2*M_PI);    

        if (ang > M_PI) {
            ang = ang - (2*M_PI);
        }
        if (ang < -M_PI) {
            ang = ang + (2*M_PI);
        }
    }
    return ang;
}

void pi_to_pi_arr(double* angles, int n) 
{
    for (int i=0; i<n; i++) {
        angles[i] = pi_to_pi(angles[i]);
    }
}
