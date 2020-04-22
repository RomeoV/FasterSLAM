#include "read_input_file.h"  // import file to test
#include <cmath>
#include <sstream>
#include <algorithm>

#include "typedefs.h"
#include "ut.hpp"
using namespace boost::ut;
using namespace boost::ut::bdd;

int main() {

"landmarks"_test = [] {
     given("I have a file with landmarks and waypoints") = [] {
        std::string filename = "example_webmap.mat";
    
        when("I read the file") = [&] {
            int num_landmarks = 35;
            int num_waypoints = 17;
            double landmarks[num_landmarks*2];
            double waypoints[num_waypoints*2];
            read_input_file(filename, landmarks, waypoints);

            then("I get arrays of landmarks and waypoints") = [&] {
                double lm_target[70] = {2.9922, -25.7009, 32.8988, -33.1776, 24.7991, -68.3801, 75.2664, -65.5763, 73.7087, -35.6698, 98.6308, 3.8941, 49.4097, 26.9470, 80.2508, 59.9688, 54.7056, 88.6293, 13.2726, 80.8411, -16.3224, 49.0654, -65.8551, 83.6449, -96.3847, 60.5919, -76.1355, 36.2928, -87.3505,-21.3396, -103.5498,-32.2430, -92.9579,-77.7259, -55.8863,-55.9190, -35.3255,-19.1589, -56.1978, 16.3551, -7.2882, 19.7819, 25.7336, 3.2710, -19.5610, 80.1444, -41.7506, 46.2325, 25.4021, 26.8543, 91.3870, 28.2385, -13.7216,-79.0338, -52.8454,-92.1833, -86.1298, 15.0890, -127.0053, 22.7018, -25.9843, 4.7078, 56.3508, -17.4388, 51.6793, -83.1863, 21.3146, -96.3358, 48.1757, 57.9979};
                double wp_target[34] = {12.6495, -41.5888, 44.7368, -54.9844, 85.5467, -45.0156, 93.6464, -17.2897, 64.9860, 5.7632, 71.8396, 31.6199, 71.2165, 70.5607, 33.5218, 76.4798, 12.9611, 51.5576, -32.5218, 67.4455, -74.8894, 69.0031, -97.3193, 41.9003, -107.5997, 6.3863, -86.4159, -25.7009, -83.9237, -64.0187, -39.3754, -81.4642, -17.5685, -51.5576};
                for(int i = 0; i<70; i++){
                    expect(that % lm_target[i] == landmarks[i]);
                }
                for(int i = 0; i<34; i++){
                    expect(that % wp_target[i] == waypoints[i]);
                }
            };
        };
    };
};

};


