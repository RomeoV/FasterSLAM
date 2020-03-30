#include "predict_true.h"

#include "pi_to_pi.h"

#include <math.h>


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

void predict_true(double V,double G,double WB,double dt, Vector3d xv) {
	xv[0] = xv[0] + V*dt*cos(G+xv[2]);		
	xv[1]= xv[1] + V*dt*sin(G+xv[2]);
	xv[2] = pi_to_pi(xv[2] + V*dt*sin(G)/WB);		
}
