#include "predict_true.h"
#include "pi_to_pi.h"
#include <math.h>

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
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
void predict_true(const double V,const double G,const double WB,
                const double dt, Vector3d xv) {
	predict_true_base(V, G, WB, dt, xv);
}
void predict_true_base(const double V,const double G,const double WB,
                const double dt, Vector3d xv) {
	xv[0] = xv[0] + V*dt*cos(G+xv[2]);		
	xv[1]= xv[1] + V*dt*sin(G+xv[2]);
	xv[2] = pi_to_pi_base(xv[2] + V*dt*sin(G)/WB);		
}
