#include <stdio.h>
#include <algorithm>
#include <iterator>
#include <unistd.h>
#include <errno.h>
#include <string>
#include <vector>

/** Reads the input file with the given path
 * Careful: This mallocs for lm and wp and expects you to dealloc!
 */
void read_input_file(const std::string s, double **lm, double **wp, size_t& N_lm, size_t& N_wp); 

/* scales by duplicating landmarks and waypoints directly, i.e.
 1 2 3 
 4 5 6 
 with scale 3 will be
 1 2 3 
 1 2 3 
 1 2 3
 4 5 6 
 4 5 6
 4 5 6 */
void read_input_file_and_scale(const std::string s, const int scale, double **lm, double **wp, size_t& N_lm, size_t& N_wp);