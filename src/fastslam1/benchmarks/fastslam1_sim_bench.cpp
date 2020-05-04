#include "rdtsc_benchmark.h"
#include "ut.hpp"

#include <random>
#include <algorithm>
#include <math.h>
#include <functional>
#include <vector>

#include "particle.h"
#include "fastslam1_utils.h"
#include "configfile.h"
#include "read_input_file.cpp"
#include "fastslam1_sim.h"


#include "predict_update.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`

int main() {
    double *lm; // landmark positions
	double *wp; // way points
    size_t lm_rows, wp_rows;

	read_input_file("../input_data/example_webmap.mat", &lm, &wp, lm_rows, wp_rows);

    Particle *particles;
	double *weights;

    Benchmark<decltype(&fastslam1_sim)> bench("fastslam1_sim Benchmark");

    bench.controls.NUM_RUNS = 1;
    bench.controls.REP = 1;
    bench.controls.CYCLES_REQUIRED = 0.0;

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&fastslam1_sim_base, "fastslam1_sim_base", 0.0);
    //bench.add_function(&fastslam1_sim_fmod, "fastslam1_sim_fmod", work);
    int N= 100;
    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    for (int i = 0; i < 5; i++) {
        NPARTICLES = pow(2,i) * N;
        bench.funcFlops[0] = NPARTICLES * 17329; // Not real work, but interesting to look at
        bench.run_name = std::to_string(NPARTICLES); // Set name of run to identify it easier
        bench.run_benchmark(lm, lm_rows, 2, wp, wp_rows, 2, &particles, &weights);
    }

    bench.details(); // We want output per func and run, so details is the choice

    cleanup_particles(&particles, &weights);

	free(lm);
	free(wp);
}
    
