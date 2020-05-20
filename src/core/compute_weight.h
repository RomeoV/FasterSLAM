#include "typedefs.h"

#include "particle.h"

double compute_weight(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]);

double compute_weight_base(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]);


double compute_weight_active(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]);

// Work / Memory instrumenting
double compute_weight_base_flops(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]);

double compute_weight_base_memory(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]);

double compute_weight_active_flops(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]);

double compute_weight_active_memory(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R,
                      Vector2d zp[],
                      Matrix23d Hv[],
                      Matrix2d Hf[],
                      Matrix2d Sf[]);
