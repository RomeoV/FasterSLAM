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

int main(int argc, char *argv[]) {
    double *lm = NULL; // landmark positions
    size_t lm_rows = 0;
    size_t N_features = 0;

    if (argc < 2){
		read_sequential_input_file("../input_data/VP_inputs.txt", &lm, lm_rows, N_features);
    } else {
        read_sequential_input_file(argv[1], &lm, lm_rows, N_features);
    }

    Particle *particles;
	double *weights;

    Benchmark<decltype(&fastslam1_sim_base_VP)> bench("fastslam1_sim_VP");

    bench.controls.NUM_RUNS = 1;
    bench.controls.REP = 1;
    bench.controls.CYCLES_REQUIRED = 0.0;

    bench.csv_path = "fastslam1_VP_particle.csv";
    bench.csv_output = false;

    // Add your functions to the struct, give it a name (Should describe improvements there) and yield the flops this function has to do (=work)
    // First function should always be the base case you want to benchmark against!
    bench.add_function(&fastslam1_sim_base_VP, "fastslam1_sim_base_VP", 0.0);
    bench.add_function(&fastslam1_sim_active_VP, "fastslam1_sim_active_VP", 0.0);
    //bench.add_function(&fastslam1_sim_fmod, "fastslam1_sim_fmod", work);
    int N= 100;
    //Run the benchmark: give the inputs of your function in the same order as they are defined. 
    for (int i = 0; i < 6; i++) {
        NPARTICLES = pow(2,i) * N;
        std::cout<< "Benchmarking N="<<NPARTICLES<<" Particles..."<<std::endl;
        bench.run_name = std::to_string(NPARTICLES); // Set name of run to identify it easier
        bench.run_benchmark(lm, lm_rows, 4, N_features, &particles, &weights);
    }

    bench.details(); // We want output per func and run, so details is the choice

    cleanup_particles_and_pose(&particles, &weights, &xv, &Pv, NPARTICLES);

    bench.write_csv_details();

	free(lm);
}