#include "typedefs.h"
#include "particle.h"

void predict_update(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_base(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_active(double *wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_simd(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
                    size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);

void predict_update_fast(double* wp, size_t N_waypoints, double V, double* Q, double dt, 
size_t N, Vector3d xtrue, int* iwp, double* G, Particle* particles);