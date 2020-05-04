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
 * Work, best: 3 sin/cos + 9 mults + 3 adds + 1 neg + 4 fl-comp +1 div = 20 flops
 * Work, worst: 3 sin/cos + 12 mults + 5 adds + 2 neg + 4 fl-comp  +2 div +1 floor = 29 flops
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
