#include <cstdlib>
#include "feature_update.h"
#include "linalg.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation
 * ToDo: Unit tests, check if correct
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: TBD
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/

void feature_update(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]) {
    feature_update_base(particle, z, idf, N_idf, R, zp, Hv, Hf, Sf);
}

// z is the list of measurements conditioned on the particle.
// void feature_update(Particle &particle, vector<Vector2d> z, vector<int>idf,
// Matrix2d R)
void feature_update_base(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]
                    ) {
  // Having selected a new pose from the proposal distribution, this pose is
  // assumed perfect and each feature update maybe computed independently and
  // without pose uncertainty

  Vector2d xf[N_idf];
  Matrix2d Pf[N_idf];

  for (size_t i = 0; i < N_idf; i++) {
    copy(particle->xf + (2 * idf[i]), 2, xf[i]);  // means
    copy(particle->Pf + (4 * idf[i]), 4, Pf[i]);  // covariances
  }

  compute_jacobians_base(particle, idf, N_idf, R, zp, Hv, Hf, Sf);

  Vector2d feat_diff[N_idf];  // difference btw feature prediciton and
                            // measurement (used to update mean)
  for (int i = 0; i < N_idf; i++) {
    sub(z[i], zp[i], 2, feat_diff[i]);
    feat_diff[i][1] = pi_to_pi_base(feat_diff[i][1]);
  }

  for (int i = 0; i < N_idf; i++) {
    KF_cholesky_update_base(xf[i], Pf[i], 
                       feat_diff[i], R, 
                       Hf[i]);
  }

  for (size_t i = 0; i < N_idf; i++) {
    set_xfi(particle, xf[i], idf[i]);
    set_Pfi(particle, Pf[i], idf[i]);
  }

}

// z is the list of measurements conditioned on the particle.
// void feature_update(Particle &particle, vector<Vector2d> z, vector<int>idf,
// Matrix2d R)
void feature_update_active(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]) {
  // Having selected a new pose from the proposal distribution, this pose is
  // assumed perfect and each feature update maybe computed independently and
  // without pose uncertainty

  Vector2d xf[N_idf] __attribute__((aligned(32)));
  Matrix2d Pf[N_idf] __attribute__((aligned(32)));

  for (size_t i = 0; i < N_idf; i++) {
    copy(particle->xf + (2 * idf[i]), 2, xf[i]);  // means
    copy(particle->Pf + (4 * idf[i]), 4, Pf[i]);  // covariances
  }

  Vector2d feat_diff[N_idf] __attribute__((aligned(32)));  // difference btw feature prediciton and
                            // measurement (used to update mean)
  for (int i = 0; i < N_idf; i++) {
    sub(z[i], zp[i], 2, feat_diff[i]);
    feat_diff[i][1] = pi_to_pi_active(feat_diff[i][1]);
  }


  for (int i = 0; i < N_idf; i++) {
    KF_cholesky_update_active(xf[i], Pf[i], 
                       feat_diff[i], R, 
                       Hf[i]);
  }

  for (size_t i = 0; i < N_idf; i++) {
    set_xfi(particle, xf[i], idf[i]);
    set_Pfi(particle, Pf[i], idf[i]);
  }
}

// Work / Memory instrumenting
double feature_update_base_flops(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]){

  double flop_count =  compute_jacobians_base(particle, idf, N_idf, R, zp, Hv, Hf, Sf) +
    N_idf * (   
      sub_flops(z[i], zp[i], 2, feat_diff[i]) +
      pi_to_pi_base_flops(feat_diff[i][1]) + 
      KF_cholesky_update_base_flops(xf[i], Pf[i], 
                       feat_diff[i], R, 
                       Hf[i])
    );
  return flop_count;
  }

double feature_update_base_memory(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]){

  double memory_called = compute_jacobians_base(particle, idf, N_idf, R, zp, Hv, Hf, Sf) + 
    N_idf * (
      copy(particle->xf + (2 * idf[i]), 2, xf[i]) +
      copy(particle->Pf + (4 * idf[i]), 4, Pf[i]) +
      sub_memory(z[i], zp[i], 2, feat_diff[i]) +
      pi_to_pi_base_memory(feat_diff[i][1]) +
      KF_cholesky_update_base(xf[i], Pf[i], 
                        feat_diff[i], R, 
                        Hf[i]) +
      set_xfi(particle, xf[i], idf[i]) +
      set_Pfi(particle, Pf[i], idf[i])
    );
  double memory_read_count = N_idf * 16;
  double memory_written_count = N_idf * 2;
  return memory_called + memory_read_count + memory_written_count;               
}

double feature_update_active_flops(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]){

  double flop_count = N_idf * (   
      sub(z[i], zp[i], 2, feat_diff[i]) +
      pi_to_pi_active_flops(feat_diff[i][1]) + 
      KF_cholesky_update_active_flops(xf[i], Pf[i], 
                       feat_diff[i], R, 
                       Hf[i])
    );
  return flop_count;
}

double feature_update_active_memory(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]){
                      
  double memory_called = N_idf * (
      copy(particle->xf + (2 * idf[i]), 2, xf[i]) +
      copy(particle->Pf + (4 * idf[i]), 4, Pf[i]) +
      sub_memory(z[i], zp[i], 2, feat_diff[i]) +
      pi_to_pi_active(feat_diff[i][1]) +
      KF_cholesky_update_active(xf[i], Pf[i], 
                        feat_diff[i], R, 
                        Hf[i]) +
      set_xfi(particle, xf[i], idf[i]) +
      set_Pfi(particle, Pf[i], idf[i])
    );
  double memory_read_count = N_idf * 16;
  double memory_written_count = N_idf * 2;
  return memory_called + memory_read_count + memory_written_count;
}
