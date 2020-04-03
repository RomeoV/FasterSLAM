#include "stratified_random.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

using namespace std;


/*****************************************************************************
 * OPTIMIZATION STATUS
 * Last Worked on: 30.03.2020
 * Done: Base Implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: TBD
 * Memory moved: TBD
 * Cycles: Not measured.
 * Performance: Not measured.
 * Optimal: Not measured.
 * Status: Not started.
 ****************************************************************************/

void stratified_random(const size_t N, double* di)
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
