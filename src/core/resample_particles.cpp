#include "stratified_resample.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include "resample_particles.h"
#include "linalg.h"

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
double find_particle_without_dependency_flops(int* count, size_t N);
double find_particle_without_dependency_memory(int* count, size_t N);
double normalize_weights_flops(double* weights, size_t N);
double normalize_weights_memory(double* weights, size_t N);
double fill_int_flops(int*x, size_t size, int val);
double fill_int_memory(int*x, size_t size, int val);
double count_occurences_flops(const size_t* indices, size_t N, int* count);
double count_occurences_memory(const size_t* indices, size_t N, int* count);


void resample_particles(Particle* particles, size_t N, double* weights,int Nmin, int doresample) {
    resample_particles_base(particles, N, weights, Nmin, doresample);
}
double resample_particles_flops(Particle* particles, size_t N, double* weights,int Nmin, int doresample) {
    return resample_particles_base_flops(particles, N, weights, Nmin, doresample);
}
double resample_particles_memory(Particle* particles, size_t N, double* weights,int Nmin, int doresample) {
    return resample_particles_base_memory(particles, N, weights, Nmin, doresample);
}

void resample_particles_base(Particle* particles, size_t N, double* weights,int Nmin, int doresample) { 
    normalize_weights(weights, N);

    double Neff;  // discard this? NO!

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
    stratified_resample_base(weights, N, &Neff, keep_indices);

    if ((Neff < Nmin) && (doresample == 1)) {
        int count[N];
        count_occurences(keep_indices,N,count);
        int i = find_particle_without_dependency(count, N);
        while (i != -1) {  // O(N^2)
            count[i] = -1;  // was 0 before
            copyParticle(particles+(keep_indices[i]), particles+i);
            count[keep_indices[i]]--;  // should make at least one particle have no dependency, so we can change it's memory

            i = find_particle_without_dependency(count, N);  // O(N)
        }
        for (int i = 0; i<N; i++) {
            weights[i] = 1.0f/N;
            particles[i].w = weights+i;
        }
    }
}
void resample_particles_base_flops_and_memory(Particle* particles, size_t N, double* weights,int Nmin, int doresample,
                                              double* FLOPS_RET, double* MEMORY_RET) { 
    double FLOPS = 0, READS = 0, WRITES = 0, MEMORY = 0;
    double* weights_internal = (double*)malloc(N * sizeof(double));
    std::copy(weights, weights+N, weights_internal);
    normalize_weights(weights_internal, N);
    FLOPS += normalize_weights_flops (weights_internal, N);
    MEMORY+= normalize_weights_memory(weights_internal, N);

    double Neff;  // discard this? NO!
    size_t keep_indices[N];  // can be seen as dependencies   
    stratified_resample_base(weights_internal, N, &Neff, keep_indices);
    FLOPS += stratified_resample_base_flops (weights_internal, N, &Neff, keep_indices);
    MEMORY+= stratified_resample_base_memory(weights_internal, N, &Neff, keep_indices);

    if ((Neff < Nmin) && (doresample == 1)) {
        int count[N];
        size_t Nfa = particles[0].Nfa;
        count_occurences(keep_indices,N,count);
        FLOPS += count_occurences_flops (keep_indices, N, count);
        MEMORY+= count_occurences_memory(keep_indices, N, count);

        int i = find_particle_without_dependency(count, N);
        FLOPS += find_particle_without_dependency_flops (count, N);
        MEMORY+= find_particle_without_dependency_memory(count, N);
        while (i != -1) {  // O(N^2)
            count[i] = -1;
            //copyParticle(particles+(keep_indices[i]), particles+i);
            FLOPS += 0; //copyParticle_flops (particles+(keep_indices[i]), particles+i);
            MEMORY+= (11 + 6 * Nfa); //copyParticle_memory(particles+(keep_indices[i]), particles+i);
            count[keep_indices[i]]--;

            i = find_particle_without_dependency(count, N);  // O(N)
            FLOPS += find_particle_without_dependency_flops (count, N);
            MEMORY+= find_particle_without_dependency_memory(count, N);
        }
        for (int i = 0; i<N; i++) {
            weights[i] = 1.0f/N;
            particles[i].w = weights+i;
        }
        FLOPS += N * (tp.div + tp.add);
        WRITES += 2;
    }
    *FLOPS_RET = FLOPS;
    *MEMORY_RET = 1 * READS + 2 * WRITES + MEMORY;
}
double resample_particles_base_flops(Particle* particles, size_t N, double* weights,int Nmin, int doresample) { 
    double flops, memory_;
    resample_particles_base_flops_and_memory(particles, N, weights, Nmin, doresample, &flops, &memory_);
    return flops;
}
double resample_particles_base_memory(Particle* particles, size_t N, double* weights,int Nmin, int doresample) { 
    double flops_, memory;
    resample_particles_base_flops_and_memory(particles, N, weights, Nmin, doresample, &flops_, &memory);
    return memory;
}

void resample_particles_orig(Particle* particles, size_t N, double* weights,int Nmin, int doresample) { 
    normalize_weights(weights, N);

    
    double Neff;  // discard this? NO!
    size_t keep_indices[N];  // can be seen as dependencies   
    stratified_resample_base(weights, N, &Neff, keep_indices);


    if ((Neff < Nmin) && (doresample == 1)) {
        int Nf = particles[0].Nf;
        Particle old_particles[N];

        for (int i = 0; i<N; i++) {
            initParticle(old_particles+i, Nf, particles[i].xv);
            copyParticle(particles+i, old_particles+i);
        }

        for (int i = 0; i< N; i++) {
            copyParticle(old_particles+keep_indices[i], particles+i);
        }
        
        for (int i = 0; i<N; i++) {
            weights[i] = 1.0f/N;
            particles[i].w = weights+i;
            delParticleMembers(old_particles+i);
        }
    }
}
double resample_particles_orig_flops(Particle* particles, size_t N, double* weights,int Nmin, int doresample) { 
    double _Neff;
    size_t* _keep;
    return normalize_weights_flops(weights, N)
           + stratified_resample_base_flops(weights, N, &_Neff, _keep)
           + N * (tp.div + tp.add);
}
double resample_particles_orig_memory(Particle* particles, size_t N, double* weights,int Nmin, int doresample) { 
    double _Neff;
    size_t* _keep;
    size_t Nfa = particles[0].Nfa;
    double copyParticle_memory = 3 * (11 + 6 * Nfa);
    return normalize_weights_memory(weights, N)
           + stratified_resample_base_memory(weights, N, &_Neff, _keep)
           + N * ( 2 * copyParticle_memory + 2 * 2 /* 2 writes */);
}

/* Basically finds a 0 */
int find_particle_without_dependency(int* count, size_t N) {
    for (int i = 0; i < N; i++) {
        if (count[i] == 0) return i;
    }
    return -1;
}
double find_particle_without_dependency_flops(int* count, size_t N) {
    int i = find_particle_without_dependency(count, N);
    if (i == -1) return N * tp.doublecomp;
    else return i * tp.doublecomp;
}
double find_particle_without_dependency_memory(int* count, size_t N) {
    int i = find_particle_without_dependency(count, N);
    if (i == -1) return N * 1;  // reads
    else return i * 1;
}


void normalize_weights(double* weights, size_t N) {
    double sum = 0;
    for (size_t i = 0; i < N; i++) {
        sum += weights[i];
    }
    for (size_t i = 0; i < N; i++) {
        weights[i] /= sum;
    }
}
double normalize_weights_flops(double* weights, size_t N) {
    return N * (tp.add + tp.div);
}
double normalize_weights_memory(double* weights, size_t N) {
    return N + 2*N;  // reads and writes
}

void fill_int(int*x, size_t size, int val) {
    assert( x != NULL );
    for (size_t i = 0; i < size; i++) {
        x[i] = val; 
    }
}
double fill_int_flops(int*x, size_t size, int val) {
    return 0;
}
double fill_int_memory(int*x, size_t size, int val) {
    return 2 * size; // writes
}

void count_occurences(const size_t* indices, size_t N, int* count) {
    fill_int(count, N, 0);

    for (size_t i = 0; i < N; i++) {
        count[indices[i]]++;
    }
}
double count_occurences_flops(const size_t* indices, size_t N, int* count) {
    return fill_int_flops(count, N, 0);
}
double count_occurences_memory(const size_t* indices, size_t N, int* count) {
    return fill_int_memory(count, N, 0) + N + 2 * N;  // func + reads + writes
}
