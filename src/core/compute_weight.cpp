#include "compute_weight.h"

#include <math.h>

#include "compute_jacobians.h"
#include "linalg.h"
#include "pi_to_pi.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation
 * ToDo: Unit tests
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

double compute_weight(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R) {
    return compute_weight_base(particle,z,N_z, idf, R);
}

double compute_weight_base(Particle* particle,
                      Vector2d z[],
                      size_t N_z,
                      int idf[],
                      Matrix2d R) {
  Vector2d zp[N_z];
  Matrix23d Hv[N_z];
  Matrix2d Hf[N_z];
  Matrix2d Sf[N_z];

  // process each feature, incrementally refine proposal distribution
  compute_jacobians_base(particle, idf, N_z, R, zp, Hv, Hf, Sf);

  Vector2d v[N_z];
  for (size_t j = 0; j < N_z; j++) {
    Vector2d v_j;
    sub(z[j], zp[j], 2, v_j);  // v_j = z[j] - zp[j]
    v_j[1] = pi_to_pi_base(v_j[1]);
    copy(v_j, 2, v[j]);  // v[j] = v_j
  }

  double w = 1.0;

  double den, num;
  // this can probably be done alot faster without this loop.....
  
  for (size_t i = 0; i < N_z; i++) {
    Matrix2d S, ST, S_inv;
    Vector2d S_inv_v;
    double vT_S_inv_v;
    // Eq. 61 in Thrun03g
    
    copy(Sf[i], 4, S);
    transpose(S, 2, 2, ST);
    inv_2x2(S, S_inv);

    
    mul(S_inv, v[i], 2, 2, 1, S_inv_v);
    mul(v[i], S_inv_v, 1, 2, 1, &vT_S_inv_v);

    

    den = 2 * M_PI * sqrt(determinant_2x2(S));
    num = exp(-0.5 * vT_S_inv_v);
    w *= (double)num / (double)den;
  }

  return w;
}
