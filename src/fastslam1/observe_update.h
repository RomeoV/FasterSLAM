#include "typedefs.h"
#include "particle.h"
#include <stddef.h>

#include <immintrin.h>

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

void observe_update_unrolled4x_plain(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) ;

void observe_update_unrolled4x(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) ;

void observe_update_simplified(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) ;

void observe_update_inplace(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights);

__m256d exp_avx2_pd (__m256d x);

void observe_update_fast_partial_avx_slower_vTMv(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights);


void observe_update_simpleavx(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights) ;

void observe_update_fast_fullavx(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights);

void observe_update_fast_KF_comp_not_unrolled(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights);

double observe_update_base_flops(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights);

double observe_update_base_memory(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights);

double observe_update_flops(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights);

double observe_update_memory(double * lm, int N_features, Vector3d xtrue, double* R, int* ftag, 
            int* da_table, int* ftag_visible, Vector2d* z, size_t* Nf_visible, Vector2d* zf, int* idf, 
            Vector2d* zn, Particle* particles, double* weights);