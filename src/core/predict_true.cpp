#include "predict_true.h"
#include "pi_to_pi.h"
#include <math.h>

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/


void predict_true(const double V,const double G,const double WB,
                const double dt, Vector3d xv) {
	predict_true_base(V, G, WB, dt, xv);
}

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work, best: 3 sin/cos + 6 mults + 3 adds + 6 pi_to_pi (est.) +1 div = 19 flops
 * Work, worst: 3 sin/cos + 6 mults + 3 adds + 15 pi_to_pi (est.) +1 div = 28 flops
 * Memory moved: 3 doubles
 * Cycles: ~ 400 cyc [170-800]
 * Performance: 0.05
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void predict_true_base(const double V,const double G,const double WB,
                const double dt, Vector3d xv) {
	xv[0] = xv[0] + V*dt*cos(G+xv[2]);		
	xv[1]= xv[1] + V*dt*sin(G+xv[2]);
	xv[2] = pi_to_pi_base(xv[2] + V*dt*sin(G)/WB);		
}
