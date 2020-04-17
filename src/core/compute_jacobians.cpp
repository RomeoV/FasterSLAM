#include "compute_jacobians.h"

#include <math.h>

#include "linalg.h"
#include "pi_to_pi.h"

void compute_jacobians(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) {
  Vector3d& xv = particle->xv;

  // vector<Vector2d> xf;
  // vector<Matrix2d> Pf;
  Vector2d xf[N_z];
  Matrix2d Pf[N_z];

  int r;
  for (size_t i = 0; i < N_z; i++) {
    copy(particle->xf + (idf[i]), 2, xf[i]);
    copy(particle->Pf + (idf[i]), 4, Pf[i]);  // particle.Pf is a array of
                                              // matrices
  }

  double dx, dy, d2, d;

  for (size_t i = 0; i < N_z; i++) {
    dx = xf[i][0] - xv[0];
    dy = xf[i][1] - xv[1];
    d2 = pow(dx, 2) + pow(dy, 2);
    d = sqrt(d2);

    Vector2d zp_vec;

    // predicted observation
    zp_vec[0] = d;
    zp_vec[1] = atan2(dy, dx) - xv[2];
    zp_vec[1] = pi_to_pi(zp_vec[1]);
    copy(zp_vec, 2, zp[i]);

    // Jacobian wrt vehicle states
    Matrix23d HvMat = {-dx / d, -dy / d, 0, dy / d2, -dx / d2, -1};

    // Jacobian wrt feature states
    Matrix2d HfMat = {dx / d, dy / d, -dy / d2, dx / d2};

    copy(HvMat, 4, Hv[4*i]);
    copy(HfMat, 4, Hf[4*i]);

    // innovation covariance of feature observation given the vehicle'
    // Eq. 60 in Thrun03g
    Matrix2d HfMat_T;
    Matrix2d Hf_Pf;
    Matrix2d Hf_Pf_HfT;
    Matrix2d Hf_Pf_HfT_R;
    transpose(HfMat, 2, 2, HfMat_T);
    mul(HfMat, Pf[i], 2, 2, 2, Hf_Pf);
    mul(Hf_Pf, HfMat_T, 2, 2, 2, Hf_Pf_HfT);
    add(Hf_Pf_HfT, R, 2 * 2, Hf_Pf_HfT_R);
    copy(Hf_Pf_HfT_R, 2 * 2, Sf[i]);
    // ........ I hate this madness
  };
}
