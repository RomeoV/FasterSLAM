#include "compute_jacobians.h"
#include <immintrin.h>

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
    //compute_jacobians_base(particle, idf, N_z, R, zp, Hv, Hf, Sf);
    compute_jacobians_fast(particle, idf, N_z, R, zp, Hv, Hf, Sf);
    //compute_jacobians_advanced_optimisations(particle, idf, N_z, R, zp, Hv, Hf, Sf);
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

    auto R_vec =  _mm256_load_pd(R);

    for (int i = 0; i < N_z; i++) {
        double px = particle->xv[0];
        double py = particle->xv[1];
        double ptheta = particle->xv[2];

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

        // Hf[i][0] = dx_dinv;
        // Hf[i][1] = dy_dinv;
        // Hf[i][2] = -dy_d2inv;
        // Hf[i][3] = dx_d2inv;

        __m256d hf_vec = _mm256_set_pd(dx_d2inv, - dy_d2inv, dy_dinv,  dx_dinv);
        _mm256_store_pd(Hf[i], hf_vec);
        //copy(HvMat, 6, Hv[i]);
        //copy(HfMat, 4, Hf[i]);
        // innovation covariance of feature observation given the vehicle'
        // Eq. 60 in Thrun03g
        // MAt x Mat
        auto pf_vec =  _mm256_load_pd(particle->Pf + 4* idf[i]);
        auto hf_perm = _mm256_permute_pd(hf_vec, 0b0101);

        auto pmm0 = _mm256_permute2f128_pd(pf_vec,pf_vec, 0b00000001); // 2 3 0 1
        auto p0p3 = _mm256_blend_pd(pf_vec, pmm0, 0b0110);
        auto p2p1 = _mm256_blend_pd(pf_vec, pmm0, 0b1001);

        auto ymm0 = _mm256_mul_pd(hf_vec, p0p3);
        auto sum = _mm256_fmadd_pd(hf_perm, p2p1, ymm0);

        //__m256d sum = _mm256_add_pd(ymm0, ymm1);

        auto hmm0 = _mm256_permute2f128_pd(hf_vec,hf_vec, 0b00000001); 
        auto h0h3 = _mm256_blend_pd(hf_vec, hmm0, 0b0110);
        auto h2h1 = _mm256_blend_pd(hf_vec, hmm0, 0b1001);

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
void compute_jacobians_nik(Particle* particle,
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


void compute_jacobians_basic_optimisations(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) {
  //std::cout << "compute_jacobians_basic_optimisations" << std::endl;
  for (size_t i = 0; i < N_z; i++) {
      //std::cout << "Nzi " << i << std::endl;

      int idf_i = idf[i];
      double *particle_xf = particle->xf;
      double *Pf_i = (particle->Pf + 4*idf_i); // length 4
      double xf_i0 = *(particle_xf + 2*idf_i);
      double xf_i1 = *(particle_xf + 2*idf_i + 1);
      
      double particle_xv_0 = particle->xv[0];
      double particle_xv_1 = particle->xv[1];
      double particle_xv_2 = particle->xv[2];

      double dx, dy, d2, d;
      dx = xf_i0 - particle_xv_0;
      dy = xf_i1 - particle_xv_1; // /xf[i][1] - particle->xv[1];
      d2 = dx * dx + dy * dy;
      d = sqrt(d2);

      // inlining of copy
      // predicted observation
      // inline/optimise atan2 and pi_to_pi_base
      double zp_vec_1 = atan2(dy, dx) - particle_xv_2;
      zp_vec_1 = pi_to_pi_base(zp_vec_1);
      zp[i][0] = d;
      zp[i][1] = zp_vec_1;

      // inlining of copy
      // Jacobian wrt vehicle states
      
      Hv[i][0] = -dx / d;
      Hv[i][1] = -dy / d;
      Hv[i][2] = 0;
      Hv[i][3] = dy / d2;
      Hv[i][4] = -dx / d2;
      Hv[i][5] = -1;

      // Jacobian wrt feature states
      // inlining of copy
      Hf[i][0] = dx / d;
      Hf[i][1] = dy / d; 
      Hf[i][2] = -dy / d2;
      Hf[i][3] = dx / d2;
      double *HfMat = Hf[i];
      
      // innovation covariance of feature observation given the vehicle'
      // Eq. 60 in Thrun03g
      Matrix2d HfMat_T;
      Matrix2d Hf_Pf;
      Matrix2d Hf_Pf_HfT;
      Matrix2d Hf_Pf_HfT_R;

      // inlining of transpose
      // transpose(HfMat, 2, 2, HfMat_T);
      // it's just a 2x2 matrix
      HfMat_T[0] = HfMat[0];
      HfMat_T[1] = HfMat[2];
      HfMat_T[2] = HfMat[1];
      HfMat_T[3] = HfMat[3];

      // inlining of mul
      // multiply 3 matrices
      // HfMat * Pf_i * HfMat_T
      // 2x2 * 2x2 * 2x2
      // so there is not much to optimise
      // this seemed to make the biggest difference
      Hf_Pf[0] = HfMat[0] * Pf_i[0] + HfMat[1] * Pf_i[2];
      Hf_Pf[1] = HfMat[0] * Pf_i[1] + HfMat[1] * Pf_i[3];
      Hf_Pf[2] = HfMat[2] * Pf_i[0] + HfMat[3] * Pf_i[2];
      Hf_Pf[3] = HfMat[2] * Pf_i[1] + HfMat[3] * Pf_i[3];
      //mul(HfMat, Pf_i, 2, 2, 2, Hf_Pf);
      Hf_Pf_HfT[0] = Hf_Pf[0] * HfMat_T[0] + Hf_Pf[1] * HfMat_T[2];
      Hf_Pf_HfT[1] = Hf_Pf[0] * HfMat_T[1] + Hf_Pf[1] * HfMat_T[3];
      Hf_Pf_HfT[2] = Hf_Pf[2] * HfMat_T[0] + Hf_Pf[3] * HfMat_T[2];
      Hf_Pf_HfT[3] = Hf_Pf[2] * HfMat_T[1] + Hf_Pf[3] * HfMat_T[3];
      //mul(Hf_Pf, HfMat_T, 2, 2, 2, Hf_Pf_HfT);

      // inlining of add
      Hf_Pf_HfT_R[0] = Hf_Pf_HfT[0] + R[0];
      Hf_Pf_HfT_R[1] = Hf_Pf_HfT[1] + R[1];
      Hf_Pf_HfT_R[2] = Hf_Pf_HfT[2] + R[2];
      Hf_Pf_HfT_R[3] = Hf_Pf_HfT[3] + R[3];

      // inlining of copy
      Sf[i][0] = Hf_Pf_HfT_R[0];
      Sf[i][1] = Hf_Pf_HfT_R[1];
      Sf[i][2] = Hf_Pf_HfT_R[2];
      Sf[i][3] = Hf_Pf_HfT_R[3];
  };

}

// around 3.5 speedup
void compute_jacobians_advanced_optimisations(Particle* particle, 
        int idf[], 
        size_t N_z,
        Matrix2d R, 
        Vector2d* zp, //measurements (range, bearing)
        Matrix23d* Hv, // jacobians of function h (deriv of h wrt pose)
        Matrix2d* Hf, // jacobians of function h (deriv of h wrt mean)
        Matrix2d* Sf) //measurement covariance
{
  // std::cout << "compute_jacobians_advanced_optimisations" << std::endl;
  // iterate over the number of features (= N_z)
  // hard to do loop unrolling because not a power of 2
  // the loops are all independent
 

  for (size_t i = 0; i < N_z; i++) {
      //std::cout << "Nzi " << i << std::endl;
      double particle_xv_0 = particle->xv[0];
      double particle_xv_1 = particle->xv[1];
      double particle_xv_2 = particle->xv[2];
      int idf_i = idf[i];
      double *particle_xf = particle->xf;
      double *Pf_i = (particle->Pf + 4*idf_i); // length 4
      double xf_i0 = *(particle_xf + 2*idf_i);
      double xf_i1 = *(particle_xf + 2*idf_i + 1);
      
      

      double dx, dy, d2, d;
      dx = xf_i0 - particle_xv_0;
      dy = xf_i1 - particle_xv_1; // /xf[i][1] - particle->xv[1];
      d2 = dx * dx + dy * dy;
      d = sqrt(d2);

      // inlining of copy
      // predicted observation
      // inline/optimise atan2 and pi_to_pi_base
      double zp_vec_1 = atan2(dy, dx) - particle_xv_2;
      zp_vec_1 = pi_to_pi_base(zp_vec_1);
      zp[i][0] = d;
      zp[i][1] = zp_vec_1;

      // Jacobian wrt feature states
      // inlining of copy
      // is slower due to _mm256_set_pds:
      //__m256d dxdy_vec = _mm256_set_pd(dx, -dy, dy, dx);
      //__m256d dd2_vec = _mm256_set_pd(d2, d2, d, d);
      //__m256d HfMat_vec = _mm256_div_pd(dxdy_vec, dd2_vec);
      // is slower too:
      //__m256d Hf_vec = _mm256_set_pd(dx / d, dy / d, -dy / d2, dx / d2);
      //_mm256_store_pd(Hf[i], HfMat_vec);

      double dxd = dx / d;
      double dyd = dy / d; 
      double dyd2 = -dy / d2;
      double dxd2 = dx / d2;
      __m256d HfMat_vec = _mm256_set_pd(dxd2, dyd2, dyd, dxd);
      _mm256_store_pd(Hf[i], HfMat_vec);

      // inlining of copy
      // Jacobian wrt vehicle states
      Hv[i][0] = - dxd;
      Hv[i][1] = - dyd;
      Hv[i][2] = 0;
      Hv[i][3] = - dyd2;
      Hv[i][4] = - dxd2;
      Hv[i][5] = -1;
      
      __m256d HfMat_T_vec = _mm256_permute4x64_pd(HfMat_vec, 0xD8);
     
      __m256d Pf_vec = _mm256_load_pd(Pf_i);
      __m256d avec = _mm256_permute4x64_pd(HfMat_vec, 0xA0); // 2,2,0,0
      __m256d bvec = _mm256_permute4x64_pd(Pf_vec, 0x44); // 1,0,1,0
      __m256d cvec = _mm256_permute4x64_pd(HfMat_vec, 0xF5); // 3,3,1,1
      __m256d dvec = _mm256_permute4x64_pd(Pf_vec, 0xEE); // 3,2,3,2
      __m256d left_mul = _mm256_mul_pd(avec,bvec);
      __m256d right_mul = _mm256_mul_pd(cvec,dvec);
      __m256d Hf_Pf_vec = _mm256_add_pd(left_mul, right_mul);

      __m256d evec = _mm256_permute4x64_pd(Hf_Pf_vec, 0xA0); // 2,2,0,0
      __m256d fvec = _mm256_permute4x64_pd(HfMat_T_vec, 0x44); // 1,0,1,0
      __m256d gvec = _mm256_permute4x64_pd(Hf_Pf_vec, 0xF5); // 3,3,1,1
      __m256d hvec = _mm256_permute4x64_pd(HfMat_T_vec, 0xEE); // 3,2,3,2
      __m256d left_mul2 = _mm256_mul_pd(evec,fvec);
      __m256d right_mul2 = _mm256_mul_pd(gvec,hvec);
      __m256d Hf_Pf_HfT_vec = _mm256_add_pd(left_mul2, right_mul2);

      __m256d R_vec = _mm256_load_pd(R);
      __m256d Hf_Pf_HfT_R_vec = _mm256_add_pd(Hf_Pf_HfT_vec, R_vec);
      _mm256_store_pd(Sf[i], Hf_Pf_HfT_R_vec);

  };

}
