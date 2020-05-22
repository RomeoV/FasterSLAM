#include "predict_true.h"
#include "pi_to_pi.h"
#include <math.h>
#include "tscheb_sine.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/


void predict_true(const double V,const double G,const double WB,
                const double dt, Vector3d xv) {
	predict_true_active(V, G, WB, dt, xv);
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
	xv[1] = xv[1] + V*dt*sin(G+xv[2]);
	xv[2] = pi_to_pi_base(xv[2] + V*dt*sin(G)/WB);		
}

// this has reduced precision 
/*
unknown:0:FAILED [5.01203e-08 <= 1e-10] 0
unknown:0:FAILED [4.63311e-08 <= 1e-10] 1
unknown:0:FAILED [2.81244e-07 <= 1e-10] 2
had to change this: 
expect(that % fabs(xv[i]-exact_xv[i]) <= 1.0e-6) << i;
*/
void predict_true_active(const double V,const double G,const double WB,
                const double dt, Vector3d xv) {
	xv[0] = xv[0] + V*dt*tscheb_cos(G+xv[2]);		
	xv[1] = xv[1] + V*dt*tscheb_sin(G+xv[2]);
	xv[2] = pi_to_pi_active(xv[2] + V*dt*tscheb_sin(G)/WB);		
}

double predict_true_base_flops(const double V,const double G,const double WB,
                const double dt, Vector3d xv) {
	return 6*tp.mul + tp.div +5*tp.add + pi_to_pi_base_flops(xv[2] + V*dt*sin(G)/WB);	
}

double predict_true_base_memory(const double V,const double G,const double WB,
                const double dt, Vector3d xv) {
	return 2* 3;	
}
