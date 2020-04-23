#pragma once
#include "particle.h"

typedef struct {
    double xtrue[3];
    double veh[2][3];
    double V;
    double alpha; // initial steer angle, used to be called G in yglee
} Vehicle;

// Setup functions
void setup_initial_Q_R();
void setup_initial_vehicle(Vehicle* v);
void setup_initial_particles(Particle **p_, double **w_, const size_t N_features);
void setup_landmarks(int **ftag_, double **da_table_, const size_t N_features);
void setup_measurements(double **z, double **zf, double **zn, int **idf, int **ftag_visible, const size_t N_features);

// Cleanup functions
void cleanup_particles(Particle **p_, double **w_);
void cleanup_landmarks(int **ftag_, double **da_table_);
void cleanup_measurements(double **z, double **zf, double **zn, int **idf, int **ftag_visible);
