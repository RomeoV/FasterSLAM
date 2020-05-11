#include "compute_jacobians.h"

#include <math.h>

#include "linalg.h"
#include "pi_to_pi.h"
#include "trigonometry.h"
#include <immintrin.h>
/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation
 * ToDo: Unit tests, check if correct
 ****************************************************************************/

void compute_jacobians(Particle* particle, int idf[], size_t N_z, Matrix2d R,
                       Vector2d zp[], Matrix23d Hv[], Matrix2d Hf[], 
                       Matrix2d Sf[]) {

    compute_jacobians_base(particle, idf, N_z, R, zp, Hv, Hf, Sf);

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

  Vector2d xf[N_z];
  Matrix2d Pf[N_z];

  int r;
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


void compute_jacobians_simd(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) {
    Vector2d xf[N_z];
    Matrix2d Pf[N_z];

    int r;
    for (size_t i = 0; i < N_z; i++) {
        copy(particle->xf + 2*idf[i], 2, xf[i]);
        copy(particle->Pf + 4*idf[i], 4, Pf[i]);  // particle.Pf is a array of matrices
    }

    double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;

    double px = particle->xv[0];
    double py = particle->xv[1];
    double ptheta = particle->xv[2];

    auto R_vec =  _mm256_load_pd(R);

    for (size_t i = 0; i < N_z; i++) {
        dx = xf[i][0] - px;
        dy = xf[i][1] - py;

        //d2 = pow(dx, 2) + pow(dy, 2);
        d2 = dx * dx + dy * dy;

        d2inv = 1.0/d2;
        d = sqrt(d2);
        dinv = 1.0/ d;

        dx_dinv = dx * dinv;
        dy_dinv = dy * dinv;
        dx_d2inv = dx * d2inv;
        dy_d2inv = dy * d2inv;

        // predicted observation
        zp[i][0] = d;
        zp[i][1] = (double)atan2(dy, dx) - ptheta;
        zp[i][1] = pi_to_pi_base(zp[i][1]);

        // Jacobian wrt vehicle states
        Matrix23d HvMat = {-dx_dinv, -dy_dinv, 0, dy_d2inv, -dx_d2inv, -1};

        // Jacobian wrt feature states
        Matrix2d HfMat = {dx_dinv, dy_dinv, -dy_d2inv, dx_d2inv};

        copy(HvMat, 6, Hv[i]);
        copy(HfMat, 4, Hf[i]);
        // innovation covariance of feature observation given the vehicle'
        // Eq. 60 in Thrun03g
        Matrix2d Hf_Pf;
        Matrix2d Hf_Pf_HfT;

        auto Pfi = *(Pf +i);

        auto pf_vec =  _mm256_load_pd(Pfi);
        auto hf_vec = _mm256_load_pd(HfMat);

        auto hf_perm = _mm256_permute_pd(hf_vec, 0b0101);

        auto pmm0 = _mm256_permute2f128_pd(pf_vec,pf_vec, 0b00000001); // 2 3 0 1
        auto hmm0 = _mm256_permute2f128_pd(hf_vec,hf_vec, 0b00000001); 

        auto p0p3 = _mm256_blend_pd(pf_vec, pmm0, 0b0110);
        auto p2p1 = _mm256_blend_pd(pf_vec, pmm0, 0b1001);

        auto h0h3 = _mm256_blend_pd(hf_vec, hmm0, 0b0110);
        auto h2h1 = _mm256_blend_pd(hf_vec, hmm0, 0b1001);

        auto ymm0 = _mm256_mul_pd(hf_vec, p0p3);
        auto ymm1 = _mm256_mul_pd(hf_perm, p2p1);

        __m256d sum = _mm256_add_pd(ymm0, ymm1);

        ymm0 = _mm256_mul_pd(h2h1, sum);
        ymm1 = _mm256_fmadd_pd(h0h3, sum, R_vec);
        ymm0 = _mm256_permute_pd(ymm0, 0b0101);
        sum = _mm256_add_pd(ymm0, ymm1);
        //print256d(sum);
        _mm256_store_pd(*(Sf +i), sum);
        //print(*Sf, 10, 2);
    };
}

void compute_jacobians_fast(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) {

    double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;

    double px = particle->xv[0];
    double py = particle->xv[1];
    double ptheta = particle->xv[2];

    auto R_vec =  _mm256_load_pd(R);

    for (int i = 0; i < N_z; i++) {
        dx = particle->xf[2*idf[i]] - px;
        dy = particle->xf[2*idf[i]+1] - py;
        //d2 = pow(dx, 2) + pow(dy, 2);
        d2 = dx * dx + dy * dy;

        
        d = sqrt(d2);
        dinv = 1.0/ d;
        d2inv = dinv * dinv;


        dx_dinv = dx * dinv;
        dy_dinv = dy * dinv;
        dx_d2inv = dx * d2inv;
        dy_d2inv = dy * d2inv;

        // predicted observation
        zp[i][0] = d;
        double theta = atan2(dy, dx) - ptheta;
        zp[i][1] = pi_to_pi(theta);

        // Jacobian wrt vehicle states
        //Matrix23d HvMat = {-dx_dinv, -dy_dinv, 0, dy_d2inv, -dx_d2inv, -1};

        Hv[i][0] = -dx_dinv;
        Hv[i][1] = -dy_dinv;
        Hv[i][2] = 0.0;
        Hv[i][3] = dy_d2inv;
        Hv[i][4] = -dx_d2inv;
        Hv[i][5] = -1.0;

        Hf[i][0] = dx_dinv;
        Hf[i][1] = dy_dinv;
        Hf[i][2] = -dy_d2inv;
        Hf[i][3] = dx_d2inv;
        //copy(HvMat, 6, Hv[i]);
        //copy(HfMat, 4, Hf[i]);
        // innovation covariance of feature observation given the vehicle'
        // Eq. 60 in Thrun03g
        // MAt x Mat
        auto pf_vec =  _mm256_load_pd(particle->Pf + 4* idf[i]);
        auto hf_vec = _mm256_load_pd(*(Hf+i));
        auto hf_perm = _mm256_permute_pd(hf_vec, 0b0101);

        auto pmm0 = _mm256_permute2f128_pd(pf_vec,pf_vec, 0b00000001); // 2 3 0 1
        auto hmm0 = _mm256_permute2f128_pd(hf_vec,hf_vec, 0b00000001); 
        auto h0h3 = _mm256_blend_pd(hf_vec, hmm0, 0b0110);
        auto h2h1 = _mm256_blend_pd(hf_vec, hmm0, 0b1001);

        auto p0p3 = _mm256_blend_pd(pf_vec, pmm0, 0b0110);
        auto p2p1 = _mm256_blend_pd(pf_vec, pmm0, 0b1001);

        auto ymm0 = _mm256_mul_pd(hf_vec, p0p3);
        auto sum = _mm256_fmadd_pd(hf_perm, p2p1, ymm0);

        //__m256d sum = _mm256_add_pd(ymm0, ymm1);

        ymm0 = _mm256_mul_pd(h2h1, sum);
        auto ymm1 = _mm256_fmadd_pd(h0h3, sum, R_vec);
        ymm0 = _mm256_permute_pd(ymm0, 0b0101);
        sum = _mm256_add_pd(ymm0, ymm1);
        //print256d(sum);
        _mm256_store_pd(*(Sf +i), sum);
        //print(*Sf, 10, 2);
    };
}

void compute_jacobians_linalg_inplace(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) {
    Vector2d xf[N_z];
    Matrix2d Pf[N_z];

    int r;
    for (size_t i = 0; i < N_z; i++) {
        copy(particle->xf + 2*idf[i], 2, xf[i]);
        copy(particle->Pf + 4*idf[i], 4, Pf[i]);  // particle.Pf is a array of matrices
    }

    double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;

    double px = particle->xv[0];
    double py = particle->xv[1];
    double ptheta = particle->xv[2];

    for (size_t i = 0; i < N_z; i++) {
        dx = xf[i][0] - px;
        dy = xf[i][1] - py;

        //d2 = pow(dx, 2) + pow(dy, 2);
        d2 = dx * dx + dy * dy;

        d2inv = 1.0/d2;
        d = sqrt(d2);
        dinv = 1.0/ d;

        dx_dinv = dx * dinv;
        dy_dinv = dy * dinv;
        dx_d2inv = dx * d2inv;
        dy_d2inv = dy * d2inv;

        // predicted observation
        zp[i][0] = d;
        zp[i][1] = atan2(dy, dx) - ptheta;
        zp[i][1] = pi_to_pi_base(zp[i][1]);

        // Jacobian wrt vehicle states
        Matrix23d HvMat = {-dx_dinv, -dy_dinv, 0, dy_d2inv, -dx_d2inv, -1};

        // Jacobian wrt feature states
        Matrix2d HfMat = {dx_dinv, dy_dinv, -dy_d2inv, dx_d2inv};

        copy(HvMat, 6, Hv[i]);
        copy(HfMat, 4, Hf[i]);
        // innovation covariance of feature observation given the vehicle'
        // Eq. 60 in Thrun03g
        Matrix2d Hf_Pf;
        Matrix2d Hf_Pf_HfT;

        Matrix2d Pfi;
        copy(*(Pf+i),4,Pfi);

        // We know Pf = symmetric!
        // Maybe also diagonal=0??


        Hf_Pf[0] = HfMat[0] * Pfi[0] +HfMat[1] * Pfi[2];
        Hf_Pf[1] = HfMat[0] * Pfi[1] +HfMat[1] * Pfi[3];
        Hf_Pf[2] = HfMat[2] * Pfi[0] +HfMat[3] * Pfi[2];
        Hf_Pf[3] = HfMat[2] * Pfi[1] +HfMat[3] * Pfi[3];

        Sf[i][0] = Hf_Pf[0] * HfMat[0] +Hf_Pf[1] * HfMat[1] ;
        Sf[i][1] = Hf_Pf[0] * HfMat[2] +Hf_Pf[1] * HfMat[3] ;
        Sf[i][2] = Hf_Pf[2] * HfMat[0] +Hf_Pf[3] * HfMat[1] ;
        Sf[i][3] = Hf_Pf[2] * HfMat[2] +Hf_Pf[3] * HfMat[3] ;

        Sf[i][0] +=R[0];
        Sf[i][1] +=R[1];
        Sf[i][2] +=R[2];
        Sf[i][3] +=R[3];
    };
}

void compute_jacobians_scalar_replacement(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) {
    Vector2d xf[N_z];
    Matrix2d Pf[N_z];

    int r;
    for (size_t i = 0; i < N_z; i++) {
        copy(particle->xf + 2*idf[i], 2, xf[i]);
        copy(particle->Pf + 4*idf[i], 4, Pf[i]);  // particle.Pf is a array of matrices
    }

    double dx, dy, d2, d, dinv, d2inv, dx_d2inv, dy_d2inv, dx_dinv, dy_dinv;

    double px = particle->xv[0];
    double py = particle->xv[1];
    double ptheta = particle->xv[2];

    for (size_t i = 0; i < N_z; i++) {
        dx = xf[i][0] - px;
        dy = xf[i][1] - py;

        //d2 = pow(dx, 2) + pow(dy, 2);
        d2 = dx * dx + dy * dy;

        d2inv = 1.0/d2;
        d = sqrt(d2);
        dinv = 1.0/ d;

        dx_dinv = dx * dinv;
        dy_dinv = dy * dinv;
        dx_d2inv = dx * d2inv;
        dy_d2inv = dy * d2inv;

        // predicted observation
        zp[i][0] = d;
        zp[i][1] = atan2(dy, dx) - ptheta;
        zp[i][1] = pi_to_pi_base(zp[i][1]);

        // Jacobian wrt vehicle states
        Matrix23d HvMat = {-dx_dinv, -dy_dinv, 0, dy_d2inv, -dx_d2inv, -1};

        // Jacobian wrt feature states
        Matrix2d HfMat = {dx_dinv, dy_dinv, -dy_d2inv, dx_d2inv};

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