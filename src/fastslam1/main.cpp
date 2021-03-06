#include "fastslam1_sim.h"
#include "read_input_file.h"
#include "particle.h"
#include "fastslam1_utils.h"
#include "configfile.h"
#include <cstdlib>


int main (int argc, char *argv[])
{

    double *lm; // landmark positions
    double *wp; // way points
    size_t lm_rows, wp_rows;

    if (argc < 2)
        return -1;

    

    Particle *particles;
    double *weights;

    int flag = 1;
    if ( argc >2 ) {
        flag = atoi(argv[2]);
    }

    if ( argc >3 ) {
        NPARTICLES = atoi(argv[3]);
    }

    read_input_file(argv[1], &lm, &wp, lm_rows, wp_rows);

    if ( flag ) {
        //Active routine
        std::cout << "Running Active" << std::endl;

        fastslam1_sim(lm, lm_rows, 2, wp, wp_rows, 2, &particles, &weights);

        cleanup_particles_and_pose(&particles, &weights, &xv, &Pv, NPARTICLES);

    } else {
        // Base routine
        std::cout << "Running Base" << std::endl;

        fastslam1_sim_base(lm, lm_rows, 2, wp, wp_rows, 2, &particles, &weights);

        cleanup_particles(&particles, &weights);
    }
    free(lm);
    free(wp);
}
