#include "read_input_file.h"  // import file to test

#include "ut.hpp"
#include <string>
#include <cmath>
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {
    //using std::filesystem::path;
    "read input file"_test = [] {
        when("I open and parse an input file with specific features") = [] {
            std::string input_file_path = "inputfiles_test/test_webmap.mat";

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
            free(lm);
            free(wp);
        };
    };

    "read input file scale"_test = [] {
        for(int scale=1; scale<6; scale++){
            when("I open and parse an input file with specific features") = [=] {
                std::string input_file_path = "inputfiles_test/test_webmap.mat";

                double *lm, *wp;
                size_t N_lm, N_wp;
                // read_input_file(input_file_path.string(), lm, wp, N_lm, N_wp);
                read_input_file_and_scale("inputfiles_test/test_webmap.mat", scale, &lm, &wp, N_lm, N_wp);
                
                then("I expect the correct number of waypoints and landmarks") = [N_lm, N_wp, scale] {
                    expect(that % N_lm == 2*scale);
                    expect(that % N_wp == 4*scale);
                };

                auto is_close = [](double lhs, double rhs) { return fabs(lhs-rhs) < 1e-14; };

                then("I expect the correct landmarks") = [lm, is_close, scale] {
                    "landmarks equality"_test = [&] () {
                        int col_lm = 3; // only for this test file
                        for (int c=0; c<scale; c++){
                            for(int i=0; i<3; i++){
                                expect(is_close(lm[c*col_lm+i], i)) << lm[c*col_lm+i] << " vs " << i << " ind " << c*col_lm+i;
                            }
                            for(int i=3; i<6; i++){
                                expect(is_close(lm[c*col_lm+i+(scale-1)*col_lm], i)) << lm[c*col_lm+i+(scale-1)*col_lm] << " vs " << i << " ind " << c*col_lm+i+(scale-1)*col_lm;
                            }
                        }
                    };
                };

                then("I expect the correct waypoints") = [wp, is_close, scale] {
                    "waypoint equality"_test = [&] () {
                        int col_wp = 2; // usually is 2
                        for (int c=0; c<scale; c++){
                           for(int i=0; i<2; i++){
                                expect(is_close(wp[c*col_wp+i], i*10)) << wp[c*col_wp+i] << " vs " << i*10 << " ind " << c*col_wp+i;
                            }
                            for(int i=2; i<4; i++){
                                expect(is_close(wp[c*col_wp+i+(scale-1)*col_wp], i*10)) << wp[c*col_wp+i+(scale-1)*3] << " vs " << i*10 << " ind " << c*col_wp+i+(scale-1)*col_wp;
                            }
                        }
                    };
                };

                free(lm);
                free(wp);
            };
        }     
    };
}