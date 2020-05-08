#include <cstdlib>

/** Tschebyscheff approximation with floats
 * Taken from here: http://mooooo.ooo/chebyshev-sine-approximation/
 */
float tscheb_fsine(float f, bool angle_is_normalized);

/** Tschebyscheff approximation with doubles
 * Taken from here: http://mooooo.ooo/chebyshev-sine-approximation/
 */
double tscheb_dsine(double f, bool angle_is_normalized);

void tscheb_dsines(double* alphas, size_t N, double* results);