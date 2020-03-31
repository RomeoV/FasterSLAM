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
    double temp = k/2;

    for (int i = 0; i<N; i++) {
        di[i] = temp + k*i + unifRand() * k - (k/2);
    }
    /*
    while (temp < 1-k/2) {
        i++;
        temp = temp+k;
        di[i] = temp + unifRand() * k - (k/2); //For-Loop below done here
    }

    //dither within interval
    vector<double>::iterator diter; 
    for (diter = di.begin(); diter != di.end(); diter++) {
        *diter = (*diter) + unifRand() * k - (k/2);
    }
    */
}

double unifRand() {
    
    return double(rand()) / double(RAND_MAX);
}
