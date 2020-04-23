#include "stratified_resample.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include "resample_particles.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base Implementation
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

void normalize_weights(double* weights, size_t N);
int find_particle_without_dependency(int* count, size_t N);
void fill_int(int *x, size_t size, int val);
void count_occurences(const size_t* indices, size_t N, int* count);

void resample_particles(Particle* particles, size_t N, double* weights)
{ 
    normalize_weights(weights, N);

    double Neff;  // discard this?

    // What is happening here:
    // Basically, we want to prevent deleting and reallocating memory,
    // since the total memory requirements after resampling are the same as before.
    // In order to resample in-place (i.e. without new allocations), we have to be careful
    // though to not override a particle with a new value while it's contents are still
    // needed for copying later.
    // 
    // Idea:
    // We basically create a DAG of which particles depend on which other particles,
    // i.e. will copy from the other particle.
    // Then we check the DAG for particles without dependencies (i.e. which will not be copied from),
    // process/fill it and remove it from the DAG, thus loosing the information of that particle.
    // We repeat this until all particles have been copied.
    // 
    // This idea relies on the fact that `keep_indices` is sorted and thus will not create
    // any circular dependencies.
    // 
    // Note that particles that kopy itself (i.e. `keep_indices[i] == i`) will stay in the
    // DAG forever and will thus not be copied --- but luckily they already contain exactly
    // what they should contain!
    size_t keep_indices[N];  // can be seen as dependencies   
    stratified_resample(weights, N, &Neff, keep_indices);

    int count[N];
    count_occurences(keep_indices, N, count);

    int i = find_particle_without_dependency(count, N);
    while (i != -1) {  // O(N^2)
        count[i] = -1;  // was 0 before
        copyParticle(&particles[keep_indices[i]], &particles[i]);
        count[keep_indices[i]]--;  // should make at least one particle have no dependency, so we can change it's memory

        i = find_particle_without_dependency(count, N);  // O(N)
    }

    /* in the end, only count[i] == 1 where keep_indices[i] = i */

    normalize_weights(weights, N);
}

/* Basically finds a 0 */
int find_particle_without_dependency(int* count, size_t N) {
    for (int i = 0; i < N; i++) {
        if (count[i] == 0) return i;
    }
    return -1;
}


void normalize_weights(double* weights, size_t N)
{
    double sum = 0;
    for (size_t i = 0; i < N; i++) {
        sum += weights[i];
    }
    for (size_t i = 0; i < N; i++) {
        weights[i] /= sum;
    }
}

void fill_int(int*x, size_t size, int val) {
    assert( x != NULL );
    for (size_t i = 0; i < size; i++) {
        x[i] = val; 
    }
}

void count_occurences(const size_t* indices, size_t N, int* count) {
    fill_int(count, N, 0);

    for (size_t i = 0; i < N; i++) {
        count[indices[i]]++;
    }
}
