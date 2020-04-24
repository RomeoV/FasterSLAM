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
                    Matrix2d R) {
    feature_update(particle, z, idf, N_idf, R);
}
// z is the list of measurements conditioned on the particle.
// void feature_update(Particle &particle, vector<Vector2d> z, vector<int>idf,
// Matrix2d R)
void feature_update_base(Particle* particle,
                    Vector2d z[],
                    int idf[],
                    size_t N_idf,
                    Matrix2d R) {
  // Having selected a new pose from the proposal distribution, this pose is
  // assumed perfect and each feature update maybe computed independently and
  // without pose uncertainty

  // double* xf = (double*)malloc(2 * N_idf * sizeof(double));
  // double* Pf = (double*)malloc(4 * N_idf * sizeof(double));
  Vector2d xf[N_idf];
  Matrix2d Pf[N_idf];

  for (size_t i = 0; i < N_idf; i++) {
    copy(particle->xf + (2 * idf[i]), 2, xf[i]);  // means
    copy(particle->Pf + (4 * idf[i]), 4, Pf[i]);  // covariances
  }

  Vector2d zp[N_idf];  // this contains all predicted features in the updated
                       // robot frame (after movement) but before considering
                       // new measurements
  Matrix23d Hv[N_idf];
  Matrix2d Hf[N_idf];
  Matrix2d Sf[N_idf];
  

  compute_jacobians_base(particle, idf, N_idf, R, zp, Hv, Hf, Sf);

  Vector2d feat_diff[N_idf];  // difference btw feature prediciton and
                            // measurement (used to update mean)
  for (int i = 0; i < N_idf; i++) {
    sub(z[i], zp[i], 2, feat_diff[i]);
    feat_diff[i][1] = pi_to_pi_base(feat_diff[i][1]);
  }

  // Vector2d vi;
  // Matrix2d Hfi;
  // Matrix2d Pfi;
  // Vector2d xfi;

  for (int i = 0; i < N_idf; i++) {
    // vi = feat_diff[i];
    // Hfi = Hf[i];
    // Pfi = Pf[i];
    // xfi = xf[i];
    // KF_cholesky_update(xfi, Pfi, vi, R, Hfi);
    KF_cholesky_update_base(xf[i], Pf[i], 
                       feat_diff[i], R, 
                       Hf[i]);
    // xf[i] = xfi;
    // Pf[i] = Pfi;
  }

  for (size_t i = 0; i < N_idf; i++) {
    set_xfi(particle, xf[i], idf[i]);
    set_Pfi(particle, Pf[i], idf[i]);
  }

  // free(feat_diff);
  // free(zp);
  // free(Hv);
  // free(Hf);
  // free(Sf);

  // free(xf);
  // free(Pf);
}
