#include "compute_steering.h"
#include "pi_to_pi.h"
#include "stdio.h"
#include <math.h>

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/
void compute_steering(cVector3d x, double* wp, const size_t N_wp, const double minD, 
                      const double rateG, const double maxG, const double dt, 
                      int* iwp, double* G) {
    compute_steering_base(x,wp,N_wp,minD,rateG,maxG,dt,iwp,G);
}


/*****************************************************************************
 * PERFORMANCE STATUS
 * Work, best: 7 p2p + 2 pow_2 + 8 add + 1 mult + 3 fl-comp + 1 atan =  22 flops
 * Work, worst: 16 p2p + 2 pow_2 + 8 add + 3 mult + 7 fl-comp + 1 atan = 37 flops 
 * 
 * #Work, best, det: 2 pow_2 + 8 add + 3 mult + 7 fl-comp + 1 atan +1 neg +1 floor =  22 flops
 * #Work, worst, det: 2 pow_2 + 10 add + 9 mult + 11 fl-comp + 1 atan +2 neg + 1 div +1 floor = 37 flops 
 * Memory moved: (3 + 2*2  +2) doubles + x ints (not counted yet)
 * Cycles: 270
 * Performance: 0.08
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void compute_steering_base(cVector3d x, double* wp, const size_t N_wp, const double minD, 
                      const double rateG, const double maxG, const double dt, 
                      int* iwp, double* G) {
    //determine if current waypoint reached
    double cwp[2];
    cwp[0] = wp[2*(*iwp)+0];
    cwp[1] = wp[2*(*iwp)+1];

    double d2 = pow((cwp[0] - x[0]),2) + pow((cwp[1]-x[1]),2);  

    if (d2 < minD*minD) {
        *iwp+=1; //switch to next
        if (*iwp >= N_wp) {
                *iwp =-1;
                return;	
        }

        cwp[0] = wp[2*(*iwp)+0];   
        cwp[1] = wp[2*(*iwp)+1];
    }


    //compute change in G to point towards current waypoint
    double deltaG = atan2(cwp[1]-x[1], cwp[0]-x[0]) - x[2] - *G;
    deltaG = pi_to_pi_base(deltaG);

    //limit rate
    double maxDelta = rateG*dt;
    if (abs(deltaG) > maxDelta) {
        int sign = (deltaG > 0) ? 1 : ((deltaG < 0) ? -1 : 0);
        deltaG = sign*maxDelta;	
    }	

    //limit angle
    double G_new = *G+deltaG;
    if (abs(G_new) > maxG) {
        int sign2 = (G_new > 0) ? 1: ((G_new < 0) ? -1 : 0);
        G_new = sign2*maxG;
    }
    *G = G_new;
}
