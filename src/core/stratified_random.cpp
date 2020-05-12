#include "stratified_random.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "fastrand.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/



void stratified_random(const size_t N, double* di) {
    stratified_random_base(N, di);
}

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: N+1 div + 2*N mult + N add (+ N random??) = 4*n+1 flops
 * Memory moved: 2*N
 * Cycles: 35*N cyc [measured with N=100]
 * Performance: 0.11
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
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

void stratified_random_fastrand(const size_t N, double* di)
{ 
    double k = 1.0/(double)N;
    //deterministic intervals

    for (int i = 0; i<N; i++) {
        di[i] = k*i + k*unifRand();
    }
}