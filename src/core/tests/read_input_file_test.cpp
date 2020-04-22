#include "read_input_file.h"  // import file to test

#include "ut.hpp"
#include <filesystem>
#include <cmath>
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    using std::filesystem::path;
    "read input file"_test = [] {
        when("I open and parse an input file with specific features") = [] {
            path input_file_path = "inputfiles_test/test_webmap.mat";

            FILE* fp = fopen("inputfiles_test/test_webmap.mat", "r"); 
            expect(fp != 0);

            double *lm, *wp;
            size_t N_lm, N_wp;
            // read_input_file(input_file_path.string(), lm, wp, N_lm, N_wp);
            read_input_file("inputfiles_test/test_webmap.mat", &lm, &wp, N_lm, N_wp);

            then("I expect the correct number of waypoints and landmarks") = [N_lm, N_wp] {
                expect(that % N_lm == 2);
                expect(that % N_wp == 4);
            };

            auto is_close = [](double lhs, double rhs) { return fabs(lhs-rhs) < 1e-14; };

            then("I expect the correct landmarks") = [lm, is_close] {
                "landmarks equality"_test = [&] (size_t i) {
                    expect(is_close(lm[i], i)) << lm[i] << " vs " << i;
                } | std::vector<size_t>{0,1,2,3,4,5};
            };

            then("I expect the correct waypoints") = [wp, is_close] {
                "waypoint equality"_test = [&] (size_t i) {
                    expect(is_close(wp[i], i*10)) << wp[i] << " vs " << i;
                } | std::vector<size_t>{0,1,2,3,4,5};
            };
        };
    };
}