#include <stdio.h>
#include <algorithm>
#include <iterator>
#include <unistd.h>
#include <errno.h>
#include <string>
#include <vector>

using namespace std;

/*!
    Stratified random vector (Sampling with subpopulations).
    @param[in]  s     path/filename, e.g. "example_webmap.mat"
    @param[out] lm    Read array of landmark positions, 35 x 2
    @param[out] wp    Read array of way points, 17 x 2
 */
void read_input_file(const string s, double *lm, double *wp); 
