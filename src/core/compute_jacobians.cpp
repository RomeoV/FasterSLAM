#include "compute_jacobians.h"

#include <math.h>

#include "linalg.h"
#include "pi_to_pi.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation
 * ToDo: Unit tests, check if correct
 ****************************************************************************/

void compute_jacobians(Particle* particle, int idf[], size_t N_z, Matrix2d R,
                       Vector2d zp[], Matrix23d Hv[], Matrix2d Hf[], 
                       Matrix2d Sf[]) {

    compute_jacobians_active(particle, idf, N_z, R, zp, Hv, Hf, Sf);

}

/*****************************************************************************
 * PERFORMANCE STATUS (N_z = number of features considered)
 * Work: N_z * (  1 add + 4 subs + 2 pow(_,2)  + 1 sqrt + 1 atan2 + 1 pi_to_pi_base \
 *              + 2 matmul 2x2x2 = 2*(8 muls + 4 adds) + 1 matadd 2x2 = 4 adds)
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void compute_jacobians_base(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) {

  Vector2d xf[N_z] __attribute__ ((aligned(32)));
  Matrix2d Pf[N_z] __attribute__ ((aligned(32)));

  for (size_t i = 0; i < N_z; i++) {
      copy(particle->xf + 2*idf[i], 2, xf[i]);
      copy(particle->Pf + 4*idf[i], 4, Pf[i]);  // particle.Pf is a array of matrices
  }

  double dx, dy, d2, d;
  
  for (size_t i = 0; i < N_z; i++) {

      dx = xf[i][0] - particle->xv[0];
      dy = xf[i][1] - particle->xv[1];
      d2 = pow(dx, 2) + pow(dy, 2);
      d = sqrt(d2);

      Vector2d zp_vec;

      // predicted observation
      zp_vec[0] = d;
      zp_vec[1] = atan2(dy, dx) - particle->xv[2];
      zp_vec[1] = pi_to_pi_base(zp_vec[1]);
      copy(zp_vec, 2, zp[i]);

      // Jacobian wrt vehicle states
      Matrix23d HvMat = {-dx / d, -dy / d, 0, dy / d2, -dx / d2, -1};

      // Jacobian wrt feature states
      Matrix2d HfMat = {dx / d, dy / d, -dy / d2, dx / d2};

      copy(HvMat, 6, Hv[i]);
      copy(HfMat, 4, Hf[i]);
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

/*****************************************************************************
 * PERFORMANCE STATUS (N_z = number of features considered)
 * Work: N_z * (  1 add + 4 subs + 2 pow(_,2)  + 1 sqrt + 1 atan2 + 1 pi_to_pi_base \
 *              + 2 matmul 2x2x2 = 2*(8 muls + 4 adds) + 1 matadd 2x2 = 4 adds)
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void compute_jacobians_active(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) {

  Vector2d xf[N_z] __attribute__ ((aligned(32)));
  Matrix2d Pf[N_z] __attribute__ ((aligned(32)));

  for (size_t i = 0; i < N_z; i++) {
      copy(particle->xf + 2*idf[i], 2, xf[i]);
      copy(particle->Pf + 4*idf[i], 4, Pf[i]);  // particle.Pf is a array of matrices
  }

  double dx, dy, d2, d;
  
  for (size_t i = 0; i < N_z; i++) {

      dx = xf[i][0] - particle->xv[0];
      dy = xf[i][1] - particle->xv[1];
      d2 = pow(dx, 2) + pow(dy, 2);
      d = sqrt(d2);

      Vector2d zp_vec;

      // predicted observation
      zp_vec[0] = d;
      zp_vec[1] = atan2(dy, dx) - particle->xv[2];
      zp_vec[1] = pi_to_pi_base(zp_vec[1]);
      copy(zp_vec, 2, zp[i]);

      // Jacobian wrt vehicle states
      Matrix23d HvMat = {-dx / d, -dy / d, 0, dy / d2, -dx / d2, -1};

      // Jacobian wrt feature states
      Matrix2d HfMat __attribute__ ((aligned(32))) = {dx / d, dy / d, -dy / d2, dx / d2};

      copy(HvMat, 6, Hv[i]);
      copy(HfMat, 4, Hf[i]);
      // innovation covariance of feature observation given the vehicle'
      // Eq. 60 in Thrun03g
      Matrix2d Hf_Pf __attribute__ ((aligned(32)));
      Matrix2d Hf_Pf_HfT __attribute__ ((aligned(32)));
      Matrix2d Hf_Pf_HfT_R;

#ifdef __AVX2__
      mm_2x2_avx_v1(HfMat, Pf[i], Hf_Pf);
      mmT_2x2_avx_v1(Hf_Pf, HfMat, Hf_Pf_HfT); 
#else      
      mm_2x2(HfMat, Pf[i], Hf_Pf);
      mmT_2x2(Hf_Pf, HfMat, Hf_Pf_HfT); 
#endif

      add(Hf_Pf_HfT, R, 2 * 2, Hf_Pf_HfT_R);
      copy(Hf_Pf_HfT_R, 2 * 2, Sf[i]);
      // ........ I hate this madness
  };

}
