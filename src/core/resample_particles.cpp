#include "stratified_resample.h"
#include <cassert>
#include <cmath>

#include "resample_particles.h"

void normalize_weights(double* weights, size_t N);
int find_particle_without_dependency(int* count, size_t N);
void fill_int(int *x, size_t size, int val);
void count_occurences(const size_t* indices, size_t N, int* count);

void resample_particles(Particle* particles, double* weights, size_t N)
{ 
    normalize_weights(weights, N);

    double Neff;  // discard this?
    size_t keep_indices[N];  // can be seen as dependencies
    stratified_resample(weights, N, Neff, keep_indices);

    int count[N];
    count_occurences(keep_indices, N, count);

    int i = find_particle_without_dependency(count, N);
    while (i != -1) {  // O(N^2)
        count[i] = -1;  // was 0 before
        copyParticle(particles[keep_indices[i]], particles[i]);
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
