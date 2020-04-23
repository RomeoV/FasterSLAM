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

void setup_initial_vehicle(Vehicle* v) {
    Vector3d xtrue_initial = {0.0, 0.0, 0.0};
    copy(v->xtrue, 3, xtrue_initial);

    double veh_initial[2][3] = { {0,-WHEELBASE,-WHEELBASE}, {0,-1,1} };
    copy(*(v->veh), 2*3, *veh_initial);

    v->V = 0; // \todo
    v->alpha = 0.0;  // along global coordinate x axis
}

void setup_initial_particles(Particle **p_, double **w_, const size_t N_features) {
    *p_ = (Particle *) malloc( NPARTICLES * sizeof(Particle) );
    Particle *p = *p_;
    for (size_t i = 0; i < NPARTICLES; i++) {
        initParticle(&p[i], N_features);
    }
    const double uniform_w = 1.0 / NPARTICLES; 
    *w_ = (double *) malloc( NPARTICLES * sizeof(double) );
    double *w = *w_;
    for (size_t i = 0; i < NPARTICLES; i++) {
        w[i] = uniform_w;
        p[i].w = &w[i];
    }
}

void cleanup_particles(Particle** p_, double** w_) {
    Particle* p = *p_;
    for (size_t i = 0; i < NPARTICLES; i++) {
        delParticleMembers(&p[i]);
    }
    free(p_);
    free(w_);
}

void setup_landmarks(int **ftag_, double **da_table_, const size_t N_features) {
    *ftag_ = (int*) malloc( N_features * sizeof(int) ); // identifier for each lm
    *da_table_ = (double *) malloc( N_features * sizeof(double) ); // data association table
    for (size_t i = 0; i < N_features; i++) {
        (*ftag_)[i] = i; 
        (*da_table_)[i] = -1;
    }
}

void cleanup_landmarks(int **ftag_, double **da_table_) {
    free(*ftag_);
    free(*da_table_);
}

// Measurements
void setup_measurements(double **z, double **zf, double **zn, int **idf, 
        int **ftag_visible, const size_t N_features) {
    *z = (double*) malloc( 2 * N_features * sizeof(double) );
    *zf = (double*) malloc( 2 * N_features * sizeof(double) );
    *zn = (double*) malloc( 2 * N_features * sizeof(double) );
    *idf = (int*) malloc( N_features * sizeof(int) );
    *ftag_visible = (int*) malloc( N_features * sizeof(int) );
}

void cleanup_measurements(double **z, double **zf, double **zn, int **idf, int **ftag_visible) {
    free(*z);
    free(*zf);
    free(*zn);
    free(*idf);
    free(*ftag_visible);
}
