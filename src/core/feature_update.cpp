#include <cstdlib>
#include "feature_update.h"
#include "linalg.h"
#include "KF_cholesky_update.h"

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
    feature_update_active(particle, z, idf, N_idf, R, zp, Hv, Hf, Sf);
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


  
  compute_jacobians_fast(particle, idf, N_idf, R, zp, Hv, Hf, Sf);

  Vector2d feat_diff;  // difference btw feature prediciton and
                            // measurement (used to update mean)
  double xf0, xf1, Pf0, Pf1, Pf2, Pf3, feat_diff0, feat_diff1;

  for (size_t i = 0; i < N_idf; i++) {
    Vector2d xf;
    Matrix2d Pf;

    xf[0] = (particle->xf + (2 * idf[i]))[0];
    xf[1] = (particle->xf + (2 * idf[i]))[1];
    Pf[0] = (particle->Pf + (4 * idf[i]))[0];
    Pf[1] = (particle->Pf + (4 * idf[i]))[1];
    Pf[2] = (particle->Pf + (4 * idf[i]))[2];
    Pf[3] = (particle->Pf + (4 * idf[i]))[3];
  
    feat_diff[0] = z[i][0] - zp[i][0];
    feat_diff1 = z[i][1] - zp[i][1];
    feat_diff[1] = pi_to_pi_active(feat_diff1);

    KF_cholesky_update_reduced_flops(xf, Pf, 
                       feat_diff, R, 
                       Hf[i]);

    particle->xf[2*idf[i]+0] = xf[0];
    particle->xf[2*idf[i]+1] = xf[1];
    particle->Pf[4*idf[i]] = Pf[0];
    particle->Pf[4*idf[i] + 1] = Pf[1];
    particle->Pf[4*idf[i] + 2] = Pf[2];
    particle->Pf[4*idf[i] + 3] = Pf[3];
  }

}

// Work / Memory instrumenting

// Careful: this function needs correct input to work, since
// pi_to_pi number of flops depends on the output of compute_jacobians
// which needs to be passed to this function via zp!
double feature_update_base_flops(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R,
                    Vector2d zp[],
                    Matrix23d Hv[],
                    Matrix2d Hf[],
                    Matrix2d Sf[]){

  Vector2d feat_diff[N_idf];
  Vector2d xf[N_idf];
  Matrix2d Pf[N_idf];

  double flop_count = compute_jacobians_base_flops(particle, idf, N_idf, R, zp, Hv, Hf, Sf) +
    N_idf * (   
      sub_flops(z[0], zp[0], 2, feat_diff[0]) +
      KF_cholesky_update_base_flops(xf[0], Pf[0], 
                       feat_diff[0], R, 
                       Hf[0])
    );

  compute_jacobians_base(particle, idf, N_idf, R, zp, Hv, Hf, Sf);
  for (int i = 0; i < N_idf; i++) {
    sub(z[i], zp[i], 2, feat_diff[i]);
    flop_count += pi_to_pi_base_flops(feat_diff[i][1]);
  }

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

  Vector2d feat_diff[N_idf];
  Vector2d xf[N_idf];
  Matrix2d Pf[N_idf];

  double memory_called = compute_jacobians_base_memory(particle, idf, N_idf, R, zp, Hv, Hf, Sf) + 
    N_idf * (
      copy_memory(particle->xf + (2 * idf[0]), 2, xf[0]) + // 2 * 2 + //copy(particle->xf + (2 * idf[i]), 2, xf[i]) +
      copy_memory(particle->Pf + (4 * idf[0]), 4, Pf[0]) + // 2 * 4 + //copy(particle->Pf + (4 * idf[i]), 4, Pf[i]) +
      sub_memory(z[0], zp[0], 2, feat_diff[0]) +
      pi_to_pi_base_memory(feat_diff[0][1]) +
      KF_cholesky_update_base_memory(xf[0], Pf[0], 
                        feat_diff[0], R, 
                        Hf[0]) +
      2 * (2 + 1) + // set_xfi(particle, xf[i], idf[i]) +
      2 * (2 + 1) // set_Pfi(particle, Pf[i], idf[i])
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

  Vector2d feat_diff[N_idf];
  Vector2d xf[N_idf];
  Matrix2d Pf[N_idf];

  double flop_count = N_idf * (   
      sub_flops(z[0], zp[0], 2, feat_diff[0]) +
      KF_cholesky_update_active_flops(xf[0], Pf[0], 
                       feat_diff[0], R, 
                       Hf[0])
    );

  for (int i = 0; i < N_idf; i++) {
    sub(z[i], zp[i], 2, feat_diff[i]);
    flop_count += pi_to_pi_active_flops(feat_diff[i][1]);
  }

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

  Vector2d feat_diff[N_idf];
  Vector2d xf[N_idf];
  Matrix2d Pf[N_idf]; 

  double memory_called = N_idf * (
      copy_memory(particle->xf + (2 * idf[0]), 2, xf[0]) + // 2 * 2 + //copy(particle->xf + (2 * idf[i]), 2, xf[i]) +
      copy_memory(particle->Pf + (4 * idf[0]), 4, Pf[0]) + //2 * 4 + //copy(particle->Pf + (4 * idf[i]), 4, Pf[i]) +
      sub_memory(z[0], zp[0], 2, feat_diff[0]) +
      pi_to_pi_active_memory(feat_diff[0][1]) +
      KF_cholesky_update_active_memory(xf[0], Pf[0], 
                        feat_diff[0], R, 
                        Hf[0]) +
      2 * (2 + 1) + // set_xfi(particle, xf[i], idf[i]) +
      2 * (2 + 1) // set_Pfi(particle, Pf[i], idf[i])
    );
  double memory_read_count = N_idf * 16;
  double memory_written_count = N_idf * 2;
  return memory_called + memory_read_count + memory_written_count;
}
