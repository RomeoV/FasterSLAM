#pragma once
#include "typedefs.h"
#include <immintrin.h>
#include "linalg.h"

/*!
    Update Kalman Filter using Cholesky decomposition. 
    @param[out] x 	    State vector.
    @param[out] P 	    Covariance matrix of state.
    @param[in]  v   	Velocity at current state.
    @param[in]  R    	Covariance matrix of measurements.
    @param[in]  H    	Jacobian of h wrt feature states.
*/
void KF_cholesky_update(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);
void KF_cholesky_update_active(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);

void KF_cholesky_update_base(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);
void KF_cholesky_update_fused_ops(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);
void KF_cholesky_update_fused_ops_avx(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);

#ifndef KF_YGLEE
void KF_cholesky_update_reduced_flops(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);
void KF_cholesky_update_reduced_flops_avx(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);

#ifdef __AVX2__
static inline void KF_cholesky_update_unrolled4_avx(__m256d *x0x2,
                                      __m256d *x1x3,
                                      __m256d *P0,
                                      __m256d *P1,
                                      __m256d *P2,
                                      __m256d *P3,
                                      __m256d const v0v2,
                                      __m256d const v1v3,
                                      __m256d const R,
                                      __m256d const H0,
                                      __m256d const H1,
                                      __m256d const H2,
                                      __m256d const H3)
{
    //! PHt = P*H^T
    __m256d const P0H0t = _mmT_2x2_avx_v3( *P0, H0 );
    __m256d const P1H1t = _mmT_2x2_avx_v3( *P1, H1 );
    __m256d const P2H2t = _mmT_2x2_avx_v3( *P2, H2 );
    __m256d const P3H3t = _mmT_2x2_avx_v3( *P3, H3 );

    //! S += H*PHt ( S = H*P*H^T + R )
    __m256d const S0 = _mmadd_2x2_avx_v2( H0, P0H0t, R );
    __m256d const S1 = _mmadd_2x2_avx_v2( H1, P1H1t, R );
    __m256d const S2 = _mmadd_2x2_avx_v2( H2, P2H2t, R );
    __m256d const S3 = _mmadd_2x2_avx_v2( H3, P3H3t, R );

    //! Sinv = S^(-1)
    __m256d S0inv, S1inv, S2inv, S3inv;
    batch_inverse_2x2(S0, S1, S2, S3, &S0inv, &S1inv, &S2inv, &S3inv);

    //! W = PHt*Sinv
    __m256d const W0 = _mm_2x2_avx_v1(P0H0t, S0inv);
    __m256d const W1 = _mm_2x2_avx_v1(P1H1t, S1inv);
    __m256d const W2 = _mm_2x2_avx_v1(P2H2t, S2inv);
    __m256d const W3 = _mm_2x2_avx_v1(P3H3t, S3inv);

//    //! Permutes for stacked mvadds
//    __m256d const x2x0 = _mm256_permute2f128_pd( *x0x2, *x0x2, 0b00000001 );
//    __m256d const v2v0 = _mm256_permute2f128_pd(  v0v2,  v0v2, 0b00000001 );
//    __m256d const x3x1 = _mm256_permute2f128_pd( *x1x3, *x1x3, 0b00000001 );
//    __m256d const v3v1 = _mm256_permute2f128_pd(  v1v3,  v1v3, 0b00000001 );
//
//    //! x = x + W*v
//    __m128d const x0 =_mvadd_2x2_avx_v1( W0, _mm256_castpd256_pd128(v0v2), _mm256_castpd256_pd128(*x0x2) );
//    __m128d const x1 =_mvadd_2x2_avx_v1( W1, _mm256_castpd256_pd128(v1v3), _mm256_castpd256_pd128(*x1x3) );
//    __m128d const x2 =_mvadd_2x2_avx_v1( W2, _mm256_castpd256_pd128(v2v0), _mm256_castpd256_pd128( x2x0) );
//    __m128d const x3 =_mvadd_2x2_avx_v1( W3, _mm256_castpd256_pd128(v3v1), _mm256_castpd256_pd128( x3x1) );
//
//    //! Store back in x0x2, x1x3
//    *x0x2 = _mm256_blend_pd( _mm256_broadcast_pd( &x0 ), _mm256_broadcast_pd( &x2 ), 0b1100 );
//    *x1x3 = _mm256_blend_pd( _mm256_broadcast_pd( &x1 ), _mm256_broadcast_pd( &x3 ), 0b1100 );

    //! x = x + W*x
    *x0x2 = _mmTadd_2x2_avx_v2(v0v2, W0, *x0x2);
    *x1x3 = _mmTadd_2x2_avx_v2(v1v3, W0, *x1x3);

    //! Z = W * (PHt)^T = P*Ht*Sinv*H*Pt
    __m256d const Z0 = _mmT_2x2_avx_v3( W0, P0H0t );
    __m256d const Z1 = _mmT_2x2_avx_v3( W1, P1H1t );
    __m256d const Z2 = _mmT_2x2_avx_v3( W2, P2H2t );
    __m256d const Z3 = _mmT_2x2_avx_v3( W3, P3H3t );

    *P0 = _mm256_sub_pd( *P0, Z0 );
    *P1 = _mm256_sub_pd( *P1, Z1 );
    *P2 = _mm256_sub_pd( *P2, Z2 );
    *P3 = _mm256_sub_pd( *P3, Z3 );
}
#endif
#endif

//! Utils
double KF_cholesky_update_base_flops(Vector2d x, Matrix2d P,
        cVector2d v, cMatrix2d R, cMatrix2d H);

double KF_cholesky_update_active_flops(Vector2d x, Matrix2d P,
        cVector2d v, cMatrix2d R, cMatrix2d H);

double KF_cholesky_update_base_memory(Vector2d x, Matrix2d P,
        cVector2d v, cMatrix2d R, cMatrix2d H);

double KF_cholesky_update_active_memory(Vector2d x, Matrix2d P,
        cVector2d v, cMatrix2d R, cMatrix2d H);
