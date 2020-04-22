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

void predict(Particle *particle, double V, double G, Matrix2d Q, double dt) {
    // Turn first, move forwards after
    // \todo Noise on input (Maybe the noise should be added at anothere place in the code?)
	if (SWITCH_PREDICT_NOISE == 1) {
		Vector2d mu = {0, 0};
		Vector2d noise;
		multivariate_gauss(mu,Q,noise);	
		V += noise[0];
		G += noise[1];
	}	
    
    double theta = particle->xv[2];
    theta = pi_to_pi(theta + G);
    particle->xv[0] += V*dt*cos(theta);
    particle->xv[1] += V*dt*sin(theta);
    particle->xv[2] = theta;
}
