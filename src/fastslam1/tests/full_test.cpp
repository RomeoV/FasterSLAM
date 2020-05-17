#include <vector>  // used for test input
#include "ut.hpp"  // import functionality and namespaces from single header file

#include "fastslam1_sim.h"
#include "read_input_file.h"
#include "particle.h"
#include "fastslam1_utils.h"
#include "configfile.h"

using namespace boost::ut;  // provides `expect`, `""_test`, etc
using namespace boost::ut::bdd;  // provides `given`, `when`, `then`


int main (int argc, char *argv[])
{
    given("Some input data that we read in") = [] {
        double *lm; // landmark positions
        double *wp; // way points
        size_t lm_rows, wp_rows;

        std::string input_file_path = "inputfiles_test/test_webmap.mat";
        read_input_file(input_file_path, &lm, &wp, lm_rows, wp_rows);


        when("We run the full sim_base") = [&] {
            Particle *particles;
            double *weights;

            fastslam1_sim_base(lm, lm_rows, 2, wp, wp_rows, 2, &particles, &weights); // TODO: Return data

            cleanup_particles(&particles, &weights);
            then("We expect it to work I guess...") = [] { ; };
        };

        when("We run the full sim") = [&] {
            Particle *particles;
            double *weights;

            fastslam1_sim(lm, lm_rows, 2, wp, wp_rows, 2, &particles, &weights); // TODO: Return data

            cleanup_particles_and_pose(&particles, &weights, &xv, &Pv, NPARTICLES);
            then("We expect it to work I guess...") = [] { ; };
        };

        free(lm);
        free(wp);
    };
}
