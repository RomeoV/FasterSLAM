#include "stratified_random.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: TBD
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/

void stratified_random(const size_t N, double* di) {
    stratified_random_base(N, di);
}
void stratified_random_base(const size_t N, double* di)
{ 
    double k = 1.0/(double)N;
    //deterministic intervals

    for (int i = 0; i<N; i++) {
        di[i] = k*i + k*unifRand();
    }
}

double unifRand() {
    
    return double(rand()) / double(RAND_MAX);
}
