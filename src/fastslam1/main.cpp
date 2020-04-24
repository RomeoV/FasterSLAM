#include "fastslam1_sim.h"
#include "read_input_file.h"
#include "particle.h"
#include "fastslam1_utils.h"

int main (int argc, char *argv[])
{
	double *lm; // landmark positions
	double *wp; // way points
    size_t lm_rows, wp_rows;

	if (argc < 2)
		return -1;

	read_input_file(argv[1], &lm, &wp, lm_rows, wp_rows);

    Particle *particles;
	double *weights;
    fastslam1_sim(lm, lm_rows, 2, wp, wp_rows, 2, &particles, &weights); // TODO: Return data


    cleanup_particles(&particles, &weights);

	free(lm);
	free(wp);

}
