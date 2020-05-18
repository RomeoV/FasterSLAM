#include "fastslam1_utils.h"
#include "configfile.h"
#include "linalg.h"
#include <math.h>
void setup_initial_Q_R() {
    // Matrix2d Q, R are declared in configfile.h
    Q[0][0] = pow(sigmaV, 2);
    Q[0][1] = 0;
    Q[1][0] = 0;
    Q[1][1] = pow(sigmaG, 2); 
    
    R[0][0] = sigmaR * sigmaR;
    R[0][1] = 0;
    R[1][0] = 0;
    R[1][1] = sigmaB * sigmaB;

    if ( SWITCH_INFLATE_NOISE == 1 ) {
        scal(*Q, 4, 2.0, *Q);
        scal(*R, 4, 2.0, *R);
    }
}


void setup_initial_particles_and_pose(Particle **p_, double **w_,  double**xv_, double**Pv_, 
                        const size_t N_particles, const size_t N_features, const Vector3d xv_initial) {
    *p_ = (Particle *) malloc( N_particles * sizeof(Particle) );
    Particle *p = *p_;
    

    // Weights and pose
    *w_ = (double *) aligned_alloc(32, N_particles * sizeof(double) );
    
    *xv_ = (double *) aligned_alloc(32, 3 * N_particles * sizeof(double) );
    *Pv_ = (double *) aligned_alloc(32, 9 * N_particles * sizeof(double) );

    double *w = *w_;
    const double uniform_w = 1.0 / N_particles; 
    for (size_t i = 0; i < N_particles; i++) {
        w[i] = uniform_w;
        p[i].w = w+i;
        p[i].xv = xv + 3*i;
        p[i].Pv = Pv + 9*i;
    }

    for (size_t i = 0; i < N_particles; i++) {
        initParticle_prealloc(p+i, N_features, xv_initial);
    }
    
}

void cleanup_particles_and_pose(Particle** p_, double** w_, double** xv_, double**Pv_, const size_t N_particles) {
    Particle* particles = *p_;
    for (size_t i = 0; i < N_particles; i++) {
        delParticleMembers_prealloc(particles+i);
    }
    free(particles);
    free(*w_);
    free(*xv_);
    free(*Pv_);
}

void setup_initial_particles(Particle **p_, double **w_, const size_t N_features, const Vector3d xv_initial) {
    const int N_particles = NPARTICLES;
    *p_ = (Particle *) malloc( N_particles * sizeof(Particle) );
    Particle *p = *p_;
    for (size_t i = 0; i < N_particles; i++) {
        initParticle(&p[i], N_features, xv_initial);
    }
    const double uniform_w = 1.0 / N_particles; 
    *w_ = (double *) malloc( N_particles * sizeof(double) );
    double *w = *w_;
    for (size_t i = 0; i < N_particles; i++) {
        w[i] = uniform_w;
        p[i].w = &w[i];
    }
}

void cleanup_particles(Particle** p_, double** w_) {
    Particle* particles = *p_;
    for (size_t i = 0; i < NPARTICLES; i++) {
        delParticleMembers(particles+i);
    }
    free(particles);
    free(*w_);
}

void setup_landmarks(int **ftag_, int **da_table_, const size_t N_features) {
    *ftag_ = (int*) malloc( N_features * sizeof(int) ); // identifier for each lm
    *da_table_ = (int *) malloc( N_features * sizeof(int) ); // data association table
    for (size_t i = 0; i < N_features; i++) {
        (*ftag_)[i] = i; 
        (*da_table_)[i] = -1;
    }
}

void cleanup_landmarks(int **ftag_, int **da_table_) {
    free(*ftag_);
    free(*da_table_);
}

// Measurements
void setup_measurements(Vector2d **z, Vector2d **zf, Vector2d **zn, int **idf, 
                        int **ftag_visible, const size_t N_features) 
{
    *z = (Vector2d*) malloc( N_features * sizeof(Vector2d) );
    *zf = (Vector2d*) malloc( N_features * sizeof(Vector2d) );
    *zn = (Vector2d*) malloc( N_features * sizeof(Vector2d) );
    *idf = (int*) malloc( N_features * sizeof(int) );
    *ftag_visible = (int*) malloc( N_features * sizeof(int) );
}

void cleanup_measurements(Vector2d **z, Vector2d **zf, Vector2d **zn, int **idf, int **ftag_visible) {
    free(*z);
    free(*zf);
    free(*zn);
    free(*idf);
    free(*ftag_visible);
}
