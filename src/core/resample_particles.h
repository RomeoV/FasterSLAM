#include "particle.h"

void resample_particles(Particle* particles, size_t N, double* weights, int Nmin, int doresample);

void resample_particles_base(Particle* particles, size_t N, double* weights, int Nmin, int doresample);

void resample_particles_orig(Particle* particles, size_t N, double* weights,int Nmin, int doresample);