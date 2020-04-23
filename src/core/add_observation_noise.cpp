#include "add_observation_noise.h"
#include "linalg.h"
#include <cstdlib>
#include <cmath>

/*****************************************************************************
 * IMPLEMENTATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 *****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: TBD
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 *****************************************************************************/

//! Adds random measurement noise. We assume R is diagnoal matrix.
//! TOCHECK: vector<Vector2d> z -> Storage choice for z, ROW-WISE for now
void add_observation_noise(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise)
{
    if ( addnoise == 1 && zlen > 0 ){
        double *randM1 = (double*)malloc(zlen*sizeof(double));
        double *randM2 = (double*)malloc(zlen*sizeof(double)); 
        fill_rand(randM1, zlen, -1.0, 1.0);
        fill_rand(randM2, zlen, -1.0, 1.0);

        const double sqrtR00 = sqrt(R[0]);
        const double sqrtR11 = sqrt(R[3]);
        for (int i = 0; i < zlen; i++) {
            z[i][0] += randM1[i]*sqrtR00; // z[i][0] = z[i][0] + randM1[i]*sqrtR00;
            z[i][1] += randM2[i]*sqrtR11; // z[i][1] = z[i][1] + randM2[i]*sqrtR11;
        }
        free(randM1);
        free(randM2);
    }	
}
