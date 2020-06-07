#include "typedefs.h"
#include "particle.h"

void predict_update(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_base(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_active(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_sine(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_old(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_simd(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_fast(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_fast_normal_rand(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

// Utils
FlopCount predict_update_base_flops(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

double predict_update_base_memory(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

FlopCount predict_update_active_flops(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

double predict_update_active_memory(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

//! ----------------------- !//
//! ---- Victoria Park ---- !//
//! ----------------------- !//

//! BASE
void predict_VP(Vector3d state, double V, double G, double *Q, double WB, double dt, bool add_control_noise);

void predict_update_VP_base(double* controls, size_t N_controls, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

//! ACTIVE
//! Takes as input the Cholesky factor S of Q ( reuse )
void predict_VP_active(Vector3d state, double V, double G, double *S, double WB, double dt, bool add_control_noise);

void predict_update_VP_active(double* controls, size_t N_controls, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

//! Takes as input the Cholesky factor S of Q ( reuse )
//! Takes as input the transformed VnGn[0] ( reuse )
void predict_VP_unrolledx4_active(Vector3d state0, Vector3d state1, Vector3d state2, Vector3d state3, 
                                    double V_, double G_, double *S, double WB, double dt, bool add_control_noise);

//! Takes as input the Cholesky factor S of Q ( reuse )
//! Takes as input the transformed VnGn[0] ( reuse )
void predict_VP_unrolledx4_active_avx(Vector3d state0, Vector3d state1, Vector3d state2, Vector3d state3, 
                                    double V_, double G_, double *S, double WB, double dt, bool add_control_noise);

