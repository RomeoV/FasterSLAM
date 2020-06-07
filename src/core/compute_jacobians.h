#pragma once
#include "typedefs.h"
#include "particle.h"

#include <math.h>

#include "linalg.h"
#include "pi_to_pi.h"
#include "tscheb_sine.h"
#include "trigonometry.h"
#include "immintrin.h"

/*!
 *  Computes the jacobians given a particle state and predict observations. [Compute-Intensive, Switch to mask]
 *  @param[in]   Particle   Particle for which the jacobian should be computed.
 *  @param[in]   idf        Feature indices.
 *  @param[in]   R          Covariance matrix of observation (diagonal).
 *  @param[out]  zp         vector of predicted observation (given the new vehicle state)
 *  @param[out]  Hv         Jacobian of h wrt vehicle states
 *  @param[out]  Hf         Jacobian of h wrt feature states
 *  @param[out]  Sf         Measurement covariance of feature observation given the vehicle.
 */

void compute_jacobians(Particle* particle, int idf[], size_t N_z, Matrix2d R,
                       Vector2d zp[], Matrix23d Hv[], Matrix2d Hf[], 
                       Matrix2d Sf[]);

void compute_jacobians_base(Particle* particle, 
        int idf[], 
        size_t N_z,
        Matrix2d R, 
        Vector2d* zp, //measurements (range, bearing)
        Matrix23d* Hv, // jacobians of function h (deriv of h wrt pose)
        Matrix2d* Hf, // jacobians of function h (deriv of h wrt mean)
        Matrix2d* Sf) //measurement covariance
;

void compute_jacobians_fast(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) ;

void compute_jacobians_fast_4particles(Particle* particle[4],
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d* zp[4],
                       Matrix23d* Hv[4],
                       Matrix2d* Hf[4],
                       Matrix2d* Sf[4]);

void compute_jacobians_fast_4particles_fullavx(Particle* particle[4],
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d* zp[4],
                       Matrix23d* Hv[4],
                       Matrix2d* Hf[4],
                       Matrix2d* Sf[4]);
                       
void compute_jacobians_active(Particle* particle, 
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]) ;
                       
void compute_jacobians_basic_optimizations(Particle* particle, 
        int idf[], 
        size_t N_z,
        Matrix2d R, 
        Vector2d* zp, //measurements (range, bearing)
        Matrix23d* Hv, // jacobians of function h (deriv of h wrt pose)
        Matrix2d* Hf, // jacobians of function h (deriv of h wrt mean)
        Matrix2d* Sf) //measurement covariance
;

void compute_jacobians_advanced_optimizations(Particle* particle, 
        int idf[], 
        size_t N_z,
        Matrix2d R, 
        Vector2d* zp, //measurements (range, bearing)
        Matrix23d* Hv, // jacobians of function h (deriv of h wrt pose)
        Matrix2d* Hf, // jacobians of function h (deriv of h wrt mean)
        Matrix2d* Sf) //measurement covariance
;

FlopCount compute_jacobians_base_flops(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]);

double compute_jacobians_base_memory(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]);

void compute_jacobians_simd(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]);



void compute_jacobians_nik(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]);


// For fastest version, i.e. compute_jacobians_fast
FlopCount compute_jacobians_active_flops(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]);

double compute_jacobians_active_memory(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]);

void compute_jacobians_scalar_replacement(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]);



void compute_jacobians_linalg_inplace(Particle* particle,
                       int idf[],
                       size_t N_z,
                       Matrix2d R,
                       Vector2d zp[],
                       Matrix23d Hv[],
                       Matrix2d Hf[],
                       Matrix2d Sf[]);

inline void compute_jacobians4x_avx(__m256d px, __m256d py,__m256d ptheta, __m256d xfp0p2, __m256d xfp1p3, 
                            __m256d Pf0, __m256d Pf1, __m256d Pf2, __m256d Pf3, __m256d R_vec,
                            __m256d* zp_dist, __m256d* zp_angle,
                            __m256d* Hf0v,  __m256d* Hf1v,  __m256d* Hf2v, __m256d* Hf3v,
                            __m256d* s_zeros, __m256d* s_ones, __m256d* s_twos, __m256d* s_threes);


static inline void compute_jacobians4x_avx_inline(__m256d px, __m256d py,__m256d ptheta, __m256d xfp0p2, __m256d xfp1p3, 
                            __m256d Pf0, __m256d Pf1, __m256d Pf2, __m256d Pf3, __m256d R_vec,
                            __m256d* zp_dist, __m256d* zp_angle,
                            __m256d* Hf0v,  __m256d* Hf1v,  __m256d* Hf2v, __m256d* Hf3v,
                            __m256d* s_zeros, __m256d* s_ones, __m256d* s_twos, __m256d* s_threes) {

    __m256d ones_intr = _mm256_set1_pd(1.);
    __m256d mones_intr = _mm256_set1_pd(-1.);
    __m256d xfi_intr, yfi_intr, dx_intr, dy_intr;
    __m256d ymm0, ymm1, ymm2, ymm3, ymm4;

    xfi_intr = _mm256_unpacklo_pd(xfp0p2, xfp1p3);
    yfi_intr = _mm256_unpackhi_pd(xfp0p2, xfp1p3);

    dx_intr = _mm256_sub_pd(xfi_intr, px);
    dy_intr = _mm256_sub_pd(yfi_intr, py);

    __m256d d2_intr  = _mm256_add_pd(_mm256_mul_pd(dx_intr, dx_intr),
                                        _mm256_mul_pd(dy_intr, dy_intr));

    __m256d d_intr = _mm256_sqrt_pd(d2_intr);
    __m256d dinv_intr = _mm256_div_pd(ones_intr, d_intr);  // _mm256_rsqrt_ps(_mm256_castpd_ps(d2_intr));
    __m256d d2inv_intr = _mm256_mul_pd(dinv_intr, dinv_intr);


    auto atan2_intr_placeholder = [](__m256d ys, __m256d xs) -> __m256d {
        __m256d result;
        for (size_t i = 0; i < 4; i++) {
            result[i] = atan2(ys[i], xs[i]);
        }
        return result;
    };

    __m256d theta_intr = _mm256_sub_pd(atan2_intr_placeholder(dy_intr, dx_intr),
                                ptheta);
    theta_intr = simple_pi_to_pi_avx(theta_intr);

    *zp_dist = d_intr;
    *zp_angle = theta_intr;

    __m256d dx_dinv_intr = _mm256_mul_pd(dx_intr, dinv_intr);
    __m256d dy_dinv_intr = _mm256_mul_pd(dy_intr, dinv_intr);
    __m256d dx_d2inv_intr = _mm256_mul_pd(dx_intr, d2inv_intr);
    __m256d dy_d2inv_intr = _mm256_mul_pd(dy_intr, d2inv_intr);

    __m256d neg_dy_d2inv = _mm256_mul_pd(mones_intr, dy_d2inv_intr);

    //Convert to vertical
    ymm0 = _mm256_permute2f128_pd(dx_dinv_intr, neg_dy_d2inv, 0b00100000);
    ymm1 = _mm256_permute2f128_pd(dx_dinv_intr, neg_dy_d2inv, 0b00110001);
    ymm2 = _mm256_permute2f128_pd(dy_dinv_intr, dx_d2inv_intr, 0b00100000);
    ymm3 = _mm256_permute2f128_pd(dy_dinv_intr, dx_d2inv_intr, 0b00110001);

    *Hf0v = _mm256_unpacklo_pd(ymm0, ymm2);
    *Hf1v = _mm256_unpackhi_pd(ymm0, ymm2);
    *Hf2v = _mm256_unpacklo_pd(ymm1, ymm3);
    *Hf3v = _mm256_unpackhi_pd(ymm1, ymm3);


    __m256d hf0_perm = _mm256_permute_pd(*Hf0v, 0b0101);
    __m256d hf1_perm = _mm256_permute_pd(*Hf1v, 0b0101);
    __m256d hf2_perm = _mm256_permute_pd(*Hf2v, 0b0101);
    __m256d hf3_perm = _mm256_permute_pd(*Hf3v, 0b0101);

    __m256d pmm0 = _mm256_permute2f128_pd(Pf0, Pf0, 0b00000001); // 2 3 0 1
    __m256d pmm1 = _mm256_permute2f128_pd(Pf1, Pf1, 0b00000001); // 2 3 0 1
    __m256d pmm2 = _mm256_permute2f128_pd(Pf2, Pf2, 0b00000001); // 2 3 0 1
    __m256d pmm3 = _mm256_permute2f128_pd(Pf3, Pf3, 0b00000001); // 2 3 0 1

    __m256d part0_pf0pf3 = _mm256_blend_pd(Pf0, pmm0, 0b0110);
    __m256d part0_pf2pf1 = _mm256_blend_pd(Pf0, pmm0, 0b1001);

    __m256d part1_pf0pf3 = _mm256_blend_pd(Pf1, pmm1, 0b0110);
    __m256d part1_pf2pf1 = _mm256_blend_pd(Pf1, pmm1, 0b1001);

    __m256d part2_pf0pf3 = _mm256_blend_pd(Pf2, pmm2, 0b0110);
    __m256d part2_pf2pf1 = _mm256_blend_pd(Pf2, pmm2, 0b1001);

    __m256d part3_pf0pf3 = _mm256_blend_pd(Pf3, pmm3, 0b0110);
    __m256d part3_pf2pf1 = _mm256_blend_pd(Pf3, pmm3, 0b1001);

    ymm0 = _mm256_mul_pd(*Hf0v, part0_pf0pf3);
    ymm1 = _mm256_mul_pd(*Hf1v, part1_pf0pf3);
    ymm2 = _mm256_mul_pd(*Hf2v, part2_pf0pf3);
    ymm3 = _mm256_mul_pd(*Hf3v, part3_pf0pf3);

    //FMA
#ifdef __FMA__
    __m256d sum0 = _mm256_fmadd_pd(hf0_perm, part0_pf2pf1, ymm0);
    __m256d sum1 = _mm256_fmadd_pd(hf1_perm, part1_pf2pf1, ymm1);
    __m256d sum2 = _mm256_fmadd_pd(hf2_perm, part2_pf2pf1, ymm2);
    __m256d sum3 = _mm256_fmadd_pd(hf3_perm, part3_pf2pf1, ymm3);
#else
    __m256d sum0 = _mm256_add_pd(_mm256_mul_pd(hf0_perm, part0_pf2pf1), ymm0);
    __m256d sum1 = _mm256_add_pd(_mm256_mul_pd(hf1_perm, part1_pf2pf1), ymm1);
    __m256d sum2 = _mm256_add_pd(_mm256_mul_pd(hf2_perm, part2_pf2pf1), ymm2);
    __m256d sum3 = _mm256_add_pd(_mm256_mul_pd(hf3_perm, part3_pf2pf1), ymm3);
#endif

    __m256d hmm0 = _mm256_permute2f128_pd(*Hf0v ,*Hf0v, 0b00000001);
    __m256d hmm1 = _mm256_permute2f128_pd(*Hf1v ,*Hf1v, 0b00000001);
    __m256d hmm2 = _mm256_permute2f128_pd(*Hf2v ,*Hf2v, 0b00000001);
    __m256d hmm3 = _mm256_permute2f128_pd(*Hf3v ,*Hf3v, 0b00000001);

    __m256d Hf0_h0h3 = _mm256_blend_pd(*Hf0v, hmm0, 0b0110);
    __m256d Hf1_h0h3 = _mm256_blend_pd(*Hf1v, hmm1, 0b0110);
    __m256d Hf2_h0h3 = _mm256_blend_pd(*Hf2v, hmm2, 0b0110);
    __m256d Hf3_h0h3 = _mm256_blend_pd(*Hf3v, hmm3, 0b0110);

    __m256d Hf0_h2h1 = _mm256_blend_pd(*Hf0v, hmm0, 0b1001);
    __m256d Hf1_h2h1 = _mm256_blend_pd(*Hf1v, hmm1, 0b1001);
    __m256d Hf2_h2h1 = _mm256_blend_pd(*Hf2v, hmm2, 0b1001);
    __m256d Hf3_h2h1 = _mm256_blend_pd(*Hf3v, hmm3, 0b1001);

    ymm0 = _mm256_mul_pd(Hf0_h2h1, sum0);
    ymm1 = _mm256_mul_pd(Hf1_h2h1, sum1);
    ymm2 = _mm256_mul_pd(Hf2_h2h1, sum2);
    ymm3 = _mm256_mul_pd(Hf3_h2h1, sum3);

#ifdef __FMA__
    __m256d xmm0 = _mm256_fmadd_pd(Hf0_h0h3, sum0, R_vec);
    __m256d xmm1 = _mm256_fmadd_pd(Hf1_h0h3, sum1, R_vec);
    __m256d xmm2 = _mm256_fmadd_pd(Hf2_h0h3, sum2, R_vec);
    __m256d xmm3 = _mm256_fmadd_pd(Hf3_h0h3, sum3, R_vec);
#else
    __m256d xmm0 = _mm256_add_pd(_mm256_mul_pd(Hf0_h0h3, sum0), R_vec);
    __m256d xmm1 = _mm256_add_pd(_mm256_mul_pd(Hf1_h0h3, sum1), R_vec);
    __m256d xmm2 = _mm256_add_pd(_mm256_mul_pd(Hf2_h0h3, sum2), R_vec);
    __m256d xmm3 = _mm256_add_pd(_mm256_mul_pd(Hf3_h0h3, sum3), R_vec);
#endif

    ymm0 = _mm256_permute_pd(ymm0, 0b0101);
    ymm1 = _mm256_permute_pd(ymm1, 0b0101);
    ymm2 = _mm256_permute_pd(ymm2, 0b0101);
    ymm3 = _mm256_permute_pd(ymm3, 0b0101);

    sum0 = _mm256_add_pd(ymm0, xmm0);
    sum1 = _mm256_add_pd(ymm1, xmm1);
    sum2 = _mm256_add_pd(ymm2, xmm2);
    sum3 = _mm256_add_pd(ymm3, xmm3);

    //Convert to horizontal
    ymm0 = _mm256_permute2f128_pd(sum0, sum2, 0b00100000);
    ymm1 = _mm256_permute2f128_pd(sum0, sum2, 0b00110001);
    ymm2 = _mm256_permute2f128_pd(sum1, sum3, 0b00100000);
    ymm3 = _mm256_permute2f128_pd(sum1, sum3, 0b00110001);

    *s_zeros = _mm256_unpacklo_pd(ymm0, ymm2);
    *s_ones = _mm256_unpackhi_pd(ymm0, ymm2);
    *s_twos = _mm256_unpacklo_pd(ymm1, ymm3);
    *s_threes = _mm256_unpackhi_pd(ymm1, ymm3);
}
