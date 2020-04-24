#pragma once
#include "particle.h"

// Setup functions
void setup_initial_Q_R();
void setup_initial_particles(Particle **p_, double **w_, const size_t N_features, const Vector3d xv_initial);
void setup_landmarks(int **ftag_, int **da_table_, const size_t N_features);
void setup_measurements(Vector2d **z, Vector2d **zf, Vector2d **zn, int **idf, int **ftag_visible, const size_t N_features);

// Cleanup functions
void cleanup_particles(Particle **p_, double **w_);
void cleanup_landmarks(int **ftag_, int **da_table_);
void cleanup_measurements(Vector2d **z, Vector2d **zf, Vector2d **zn, int **idf, int **ftag_visible);
