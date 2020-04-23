#include "fastslam1_sim.h"
#include "read_input_file.h"
#include "particle.h"
#include "fastslam1_utils.h"

int main (int argc, char *argv[])
{
//	MyTimer Timer = MyTimer();
//	Timer.Start();
	double *lm; // landmark positions
	double *wp; // way points
    size_t lm_rows, wp_rows;

	if (argc < 2)
		return -1;

	read_input_file(argv[1], &lm, &wp, lm_rows, wp_rows);
//    for (int i = 0; i < lm_rows; i++) {
//        std::cout << lm[i*2+0] << " " << lm[i*2+1] << std::endl;
//    }
//    std::cout << std::endl;
//    for (int i = 0; i < wp_rows; i++) {
//        std::cout << wp[i*2+0] << " " << wp[i*2+1] << std::endl;
//    }
    Particle *particles;
	double *weights;
    fastslam1_sim(lm, lm_rows, 2, wp, wp_rows, 2, &particles, &weights); // TODO: Return data
/*
	for (int i = 0; i < data.size(); i++) {
		std::cout << "particle i=" << i << std::endl;
		std::cout << std::endl;
		std::cout << "weight" << std::endl;
		std::cout << data[i].w() << std::endl;
		std::cout << std::endl;
		std::cout << "xv (robot pose)" << std::endl;
		std::cout << data[i].xv() << std::endl;
		std::cout << std::endl;
		std::cout << "Pv (controls)" << std::endl;
		std::cout << data[i].Pv() << std::endl;
		std::cout << std::endl;
		std::cout << "xf (EFK means)" << std::endl;
		for(int j=0; j<data[i].xf().size(); j++) {
			std::cout << data[i].xf()[j] << std::endl;
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << "Pf (covariance mat)" << std::endl;
		for(int k = 0; k < data[i].Pf().size(); k++) {
			std::cout << data[i].Pf()[k] << std::endl;
			std::cout << std::endl;
		}
		std::cout << std::endl;
		std::cout << std::endl;
	}
*/
//	Timer.Stop();
//	Timer.Print("fastslam 1.0 ");
    cleanup_particles(&particles, &weights);
}
