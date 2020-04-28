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