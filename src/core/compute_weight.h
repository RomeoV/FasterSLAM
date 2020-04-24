#include "typedefs.h"

#include "particle.h"

double compute_weight(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R);

double compute_weight_base(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R);