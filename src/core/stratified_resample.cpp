#include "stratified_resample.h"
#include "stratified_random.h"
#include "linalg.h"

#include <math.h>
#include <stdio.h>

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base Implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work, best: (4*N +1) str.ran. + 3N adds + N mults + N divs + N fl.comp
 * Work, worst: (4*N +1) str.ran. + 3N adds + N mults + N divs + N^2 fl.comp
 * Memory moved: 2*N doubles + N ints + 2*N doubles
 * Cycles: 100 * N (9*n <-> 11*N)
 * Performance: 0.08
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void stratified_resample(const double* w_, const size_t N_w, double* Neff, size_t* keep) {
    stratified_resample_base(w_, N_w, Neff, keep);
}
double stratified_resample_flops(const double* w_, const size_t N_w, double* Neff, size_t* keep) {
    return stratified_resample_base_flops(w_, N_w, Neff, keep);
}
double stratified_resample_memory(const double* w_, const size_t N_w, double* Neff, size_t* keep) {
    return stratified_resample_base_flops(w_, N_w, Neff, keep);
}

void stratified_resample_base(const double* w_, const size_t N_w, double* Neff, size_t* keep) {
    // create local copy
    double *w = (double *) malloc( N_w * sizeof(double) ); 
    copy(w_, N_w, w);

    double wsum = 0.0;
    double wsqrd_sum = 0.0;

    for (int i = 0; i< N_w; i++){
        wsum += w[i];
    }

    for (int i=0; i<N_w; i++) {
        w[i] = w[i]/ wsum;
        wsqrd_sum += w[i] * w[i];
    }
    
    *Neff = 1.0f/wsqrd_sum;

    double *select = (double*) malloc(N_w * sizeof(double));

    stratified_random(N_w,select); 

    cumsum_base(w, N_w); 

    int ctr=0;
    for (int i=0; i<N_w; i++) {
        while ((ctr<N_w) && (select[ctr]<w[i])) {
            keep[ctr] = i;
            ctr++;
        }
    }
    

    free(select);
    free(w);    
}
double stratified_resample_base_flops(const double* w_, const size_t N_w, double* Neff, size_t* keep) {
    double * _;
    return N_w * ( 2*tp.add + tp.div + tp.mul + tp.doublecomp ) + tp.div + stratified_random_flops(N_w, _) + cumsum_base_flops(_, N_w);
}
double stratified_resample_base_memory(const double* w_, const size_t N_w, double* Neff, size_t* keep) {
    double *_;
    double READS = N_w * ( 5 );
    double WRITES = N_w * ( 3 /* one is integer writes */ );
    return READS + 2*WRITES + stratified_random_memory(N_w, _) + cumsum_base_memory(_, N_w);
}

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Last Worked on: 30.03.2020
 * Done: Base Implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: N adds
 * Memory moved: 2*N
 * Cycles: Not measured.
 * Performance: Not measured.
 * Optimal: Not measured.
 * Status: Not started.
 ****************************************************************************/

void cumsum(double* w, const size_t N_w) {
    cumsum_base(w, N_w);
}
double cumsum_flops(double* w, const size_t N_w) {
    return cumsum_base_flops(w, N_w);
}
double cumsum_memory(double* w, const size_t N_w) {
    return cumsum_base_memory(w, N_w);
}

void cumsum_base(double* w, const size_t N_w) {
    for (int i = 1; i<N_w; i++) {
        w[i]+=w[i-1];
    }
}
double cumsum_base_flops(double* w, const size_t N_w) {
    return (N_w - 1)*tp.add;
}
double cumsum_base_memory(double* w, const size_t N_w) {
    return N_w-1;
}

