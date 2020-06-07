#include "particle.h"

void resample_particles(Particle* particles, size_t N, double* weights, int Nmin, int doresample);
FlopCount resample_particles_flops(Particle* particles, size_t N, double* weights, int Nmin, int doresample);
double resample_particles_memory(Particle* particles, size_t N, double* weights, int Nmin, int doresample);

void resample_particles_dag(Particle* particles, size_t N, double* weights,int Nmin, int doresample);
FlopCount resample_particles_dag_flops(Particle* particles, size_t N, double* weights,int Nmin, int doresample);
double resample_particles_dag_memory(Particle* particles, size_t N, double* weights,int Nmin, int doresample);

void resample_particles_base(Particle* particles, size_t N, double* weights, int Nmin, int doresample);
FlopCount resample_particles_base_flops(Particle* particles, size_t N, double* weights, int Nmin, int doresample);
double resample_particles_base_memory(Particle* particles, size_t N, double* weights, int Nmin, int doresample);
