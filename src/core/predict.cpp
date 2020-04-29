#include "configfile.h"
#include "predict.h"
#include "pi_to_pi.h"
#include "multivariate_gauss.h"
#include <cmath>

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

void predict(Particle *particle, double V, double G, Matrix2d Q, double WB, double dt) {
	predict_base(particle, V, G, Q, WB, dt);
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
    particle->xv[0] += V*dt*cos(xv2);
    particle->xv[1] += V*dt*sin(xv2); 
    particle->xv[2] = pi_to_pi(xv2 + V*dt*sin(G)/WB);
}
