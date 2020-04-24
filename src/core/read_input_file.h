#include <stdio.h>
#include <algorithm>
#include <iterator>
#include <unistd.h>
#include <errno.h>
#include <string>
#include <vector>

using namespace std;
/*!
    Reads the input file with the given path.
    Careful: This mallocs for lm and wp and expects you to dealloc!
    @param[in]  s     path/filename, e.g. "example_webmap.mat"
    @param[out] lm    Read array of landmark positions, 35 x 2
    @param[out] wp    Read array of way points, 17 x 2
 */
void read_input_file(const string s, double **lm, double **wp, size_t& N_lm, size_t& N_wp); 