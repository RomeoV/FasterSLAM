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
 * Work: TBD
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/

void stratified_resample(const double* w_, const size_t N_w, double* Neff, size_t* keep) {
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

    cumsum(w, N_w); 

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

void cumsum(double* w, const size_t N_w) {
    for (int i = 1; i<N_w; i++) {
        w[i]+=w[i-1];
    }
}

