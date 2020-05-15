#include "typedefs.h"
#include "particle.h"
#include <stddef.h>

void observe_update(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, 
            int* idf, Vector2d* zn, Particle* particles, double* weights);

void observe_update_base(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, 
            int* idf, Vector2d* zn, Particle* particles, double* weights);

void observe_update_active(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, 
            int* idf, Vector2d* zn, Particle* particles, double* weights);

void observe_update_fast(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, 
            int* idf, Vector2d* zn, Particle* particles, double* weights);
