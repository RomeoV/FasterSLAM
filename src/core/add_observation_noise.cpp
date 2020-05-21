#include "add_observation_noise.h"
#include "linalg.h"
#include <cstdlib>
#include <cmath>

/*****************************************************************************
 * IMPLEMENTATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 *****************************************************************************/



//! Adds random measurement noise. We assume R is diagnoal matrix.
//! TOCHECK: vector<Vector2d> z -> Storage choice for z, ROW-WISE for now

void add_observation_noise(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise) {
    add_observation_noise_active(z, zlen, R, addnoise);
}
double add_observation_noise_flops(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise) {
    return add_observation_noise_base_flops(z, zlen, R, addnoise);
}
double add_observation_noise_memory(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise) {
    return add_observation_noise_base_memory(z, zlen, R, addnoise);
}


/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: 2*zlen rand + 2 sqrt + 2*zlen adds + 2*zlen mults = 6*zlen + 2
 * Memory moved: TBD
 * Cycles: 350
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 *****************************************************************************/
void add_observation_noise_base(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise)
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

void add_observation_noise_active(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise)
{
    if ( addnoise == 1 && zlen > 0 ){
        const double sqrtR00 = sqrt(R[0]);
        const double sqrtR11 = sqrt(R[3]);
        double randM1i, randM2i, randM1i_sqrtR00, randM2i_sqrtR11;
        for (int i = 0; i < zlen; i++) {
            //__m256d rand_vec1 = fill_rand_avx(-1.0,1.0);
            // inline fill_rand
            randM1i = -1.0 + 2 * ((double)rand())/((double)RAND_MAX);
            randM2i = -1.0 + 2 * ((double)rand())/((double)RAND_MAX);
            randM1i_sqrtR00 = randM1i * sqrtR00;
            randM2i_sqrtR11 = randM2i * sqrtR11;
            z[i][0] += randM1i_sqrtR00;
            z[i][1] += randM2i_sqrtR11;
        }
    }	
}


double add_observation_noise_base_flops(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise) {
    if( addnoise == 1 && zlen > 0 ) {
        double* _m1;
        return 2*fill_rand_flops(_m1, zlen, -1, 1) + 2 * tp.sqrt + zlen * 2 * ( tp.mul + tp.add );
    }
    else return 0;

}
double add_observation_noise_base_memory(Vector2d z[], const size_t zlen, cMatrix2d R, const int addnoise) {
    if( addnoise == 1 && zlen > 0 ) {
        double* _m1;
        return 2*fill_rand_memory(_m1, zlen, -1, 1) + 2 + zlen * (2 + 2*2);
    }
    else return 0;
}