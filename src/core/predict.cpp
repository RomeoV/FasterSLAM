#include "configfile.h"
#include "predict.h"
#include "pi_to_pi.h"
#include "multivariate_gauss.h"
#include <cmath>
#include <iostream>
#include "tscheb_sine.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work, best: 20 (predict_true) = 20 flop
 * Work, worst: 29 (predict_true) + 17 (multivar) = 46 flop
 * 
 * #work best, det: 3 sin/cos + 9 mults + 3 adds + 1 neg + 4 fl-comp +1 div = 20 flops
 * #work, worst, det: 3 neg + 17 mults + 4 div + 1 floor + 4 fl-comp + 3 sin + 12 adds + 2 sqrt = 46 flops
 * Memory moved: 5 doubles
 * Cycles: 230
 * Performance: 0.07
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/

void predict(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt) {
	predict_active(particle, V, G, Q, WB, dt);
}

void predict_base(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt) {
    // Turn first, move forwards after
    // \todo Noise on input (Maybe the noise should be added at anothere place in the code?)
	if (SWITCH_PREDICT_NOISE == 1) {
		Vector2d mu = {V, G};
		Vector2d noise;
		multivariate_gauss_base(mu,Q,noise);	
		V = noise[0];
		G = noise[1];
	}	
    
    double xv2 = particle->xv[2];
    particle->xv[0] += V*dt*cos(G + xv2);
    particle->xv[1] += V*dt*sin(G + xv2); 
    particle->xv[2] = pi_to_pi_base(xv2 + V*dt*sin(G)/WB);
}

// this has reduced precision 
/*
unknown:0:FAILED [3.17282e-09 <= 1e-10] 0
unknown:0:FAILED [6.05159e-09 <= 1e-10] 1
unknown:0:FAILED [5.33904e-08 <= 1e-10] 2
had to change this: 
expect(that % fabs(xv[i]-exact_xv[i]) <= 1.0e-6) << i;
*/
void predict_active(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt) {
    // Turn first, move forwards after
    // \todo Noise on input (Maybe the noise should be added at anothere place in the code?)
	if (SWITCH_PREDICT_NOISE == 1) {
		Vector2d mu = {V, G};
		Vector2d noise;
		multivariate_gauss_active(mu,Q,noise);	
		V = noise[0];
		G = noise[1];
	}	
    
    double xv2 = particle->xv[2];
    particle->xv[0] += V*dt*tscheb_cos(G + xv2);
    particle->xv[1] += V*dt*tscheb_sin(G + xv2); 
    particle->xv[2] = pi_to_pi_active(xv2 + V*dt*tscheb_sin(G)/WB);
}

double predict_base_flops(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt) {
	double flop_count = 0.0;
	if (SWITCH_PREDICT_NOISE == 1) {
		Vector2d mu = {V, G};
		Vector2d noise;
		flop_count+= multivariate_gauss_base_flops(mu, Q, noise);
	}
	flop_count+= 6*tp.mul + tp.div +5*tp.add + pi_to_pi_base_flops(particle->xv[2] + V*dt*sin(G)/WB);
	return flop_count;
}

double predict_base_memory(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt) {
	double memory_moved = 0.0;

	if (SWITCH_PREDICT_NOISE == 1) {
		Vector2d mu = {V, G};
		Vector2d noise;
		memory_moved+= multivariate_gauss_base_memory(mu, Q, noise);
	}
	return memory_moved + 2*3.0;
}

