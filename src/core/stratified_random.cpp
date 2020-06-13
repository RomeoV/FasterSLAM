#include "stratified_random.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "fastrand.h"
#include "typedefs.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/



void stratified_random(const size_t N, double* di) {
    stratified_random_base(N, di);
}
double stratified_random_flops(const size_t N, double* di) {
    return stratified_random_base_flops(N, di);
}
double stratified_random_memory(const size_t N, double* di) {
    return stratified_random_base_memory(N, di);
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
double stratified_random_base_flops(const size_t N, double* di) { 
    return tp.div + N * ( 2*tp.mul + tp.add + unifRand_flops() );
}
double stratified_random_base_memory(const size_t N, double* di) { 
    double WRITES = N;
    return 2 * WRITES;
}

double unifRand() {
    
    return double(rand()) / double(RAND_MAX);
}
double unifRand_flops() {
    return tp.rand + tp.div;
}
double unifRand_memory() {
    return 0;
}

// THIS CURRENTLY JUST CALLS REGULAR RAND!
void stratified_random_fastrand(const size_t N, double* di)
{ 
    double k = 1.0/(double)N;
    //deterministic intervals

    for (int i = 0; i<N; i++) {
        di[i] = k*i + k*unifRand();
    }
}
double stratified_random_fastrand_flops(const size_t N, double* di) {
    return tp.div + N * ( 2*tp.mul + tp.add + tp.fastrand );
}
double stratified_random_fastrand_memory(const size_t N, double* di) {
    double WRITES = N;
    return 2 * WRITES;
}