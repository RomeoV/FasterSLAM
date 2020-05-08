#include <cstdlib>

/** Tschebyscheff approximation with floats
 * Taken from here: http://mooooo.ooo/chebyshev-sine-approximation/
 */
float tscheb_fsine(float f, bool angle_is_normalized);

/** Tschebyscheff approximation with doubles
 * Taken from here: http://mooooo.ooo/chebyshev-sine-approximation/
 */
double tscheb_dsine(double f, bool angle_is_normalized);

/** Looped version of the others
 * Assumes angles are normalized! If not, you need to normalize them first (using pi_to_pi.h)
 */
void tscheb_dsines(double* alphas, size_t N, double* results);

/** Unrolled version of tscheb_dsines(double*, size_t, double*)
 * Assumes angles are normalized! If not, you need to normalize them first (using pi_to_pi.h)
 * Only a little bit faster than the non-unrolled version
 */
void tscheb_dsines_unrolled(double* alphas, size_t N, double* results);

void tscheb_dsines_avx(double* alphas, size_t N, double* results);