#include "fastslam1_sim.h"
#include "read_input_file.h"
#include "particle.h"
#include "fastslam1_utils.h"
#include "configfile.h"
#include <cstdlib>


int main (int argc, char *argv[])
{

	double *lm = NULL; // landmark positions
    size_t lm_rows = 0;
	size_t N_features = 0;

	if (argc < 2)
		return -1;

    read_sequential_input_file(argv[1], &lm, lm_rows, N_features);

	Particle *particles;
	double *weights;

    int flag = 1;
    if ( argc == 3 ) {
        flag = atoi(argv[2]);
    }

//  // deactivate for now since there is no active implementation
//	if ( flag ) {
//		//Active routine
//		fastslam1_sim_VP(lm, lm_rows, 4, N_features, &particles, &weights);
//
//		cleanup_particles_and_pose(&particles, &weights, &xv, &Pv, NPARTICLES);
//
//	} else {
		// Base routine
		fastslam1_sim_base_VP(lm, lm_rows, 4, N_features, &particles, &weights);

		cleanup_particles(&particles, &weights);
//	}
	free(lm);
}
