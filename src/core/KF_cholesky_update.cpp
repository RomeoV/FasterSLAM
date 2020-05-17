#include "KF_cholesky_update.h"
#include "linalg.h"
#include <immintrin.h>

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

void KF_cholesky_update(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H) {
    // KF_cholesky_update_base(x, P, v, R, H);
#ifdef KF_YGLEE
    #ifdef __AVX2__
        KF_cholesky_update_fused_ops_avx(x, P, v, R, H);
    #else
        KF_cholesky_update_fused_ops(x, P, v, R, H);
    #endif
#else
    #ifdef __AVX2__
        KF_cholesky_update_reduced_flops_avx(x, P, v, R, H);
    #else
        KF_cholesky_update_reduced_flops(x, P, v, R, H);
    #endif
#endif
}

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: 2 neg + 60 adds + 54 muls + 2 divs + 2 sqrts = 120 flops
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void KF_cholesky_update_base(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H)
{
    double Ht[4], PHt[4], HPHt[4];
    double S[4], St[4], SChol[4], SCholInv[4], SCholInvt[4];
    double W1[4], W[4], Wv[2], W1t[4], W1W1t[4];

    transpose(H, 2, 2, Ht);               //! Ht = H.transpose();
    mul(P, Ht, 2, 2, 2, PHt);             //! PHt = P*Ht;
    mul(H, PHt, 2, 2, 2, HPHt);           //! HPHt = H*PHt;
    add(HPHt, R, 4, S);                   //! S = HPHt + R;

    transpose(S, 2, 2, St);               //! St = S.transpose();
    add(S, St, 4, S);                     //! S = S + St
    scal(S, 4, 0.5, S);                   //! S = 0.5*S;

    llt_2x2(S, SChol);                    //! S = SChol * SChol^T
    inv_2x2(SChol, SCholInv);             //! SCholInv = inv( SChol )

#ifdef KF_YGLEE
    mul(PHt, SCholInv, 2, 2, 2, W1);     //! W1 = PHt * SCholInv;
    transpose(SCholInv, 2, 2, SCholInvt); //! SCholInvt = SCholInv.transpose();
    mul(W1, SCholInvt, 2, 2, 2, W);        //! W = W1 * SCholInvt;
#else
    transpose(SCholInv, 2, 2, SCholInvt); //! SCholInvt = SCholInv.transpose();
    mul(PHt, SCholInvt, 2, 2, 2, W1);     //! W1 = PHt * SCholInvt;
    mul(W1, SCholInv, 2, 2, 2, W);        //! W = W1 * SCholInv;
#endif

    //! x = x + Wv;
    mul(W, v, 2, 2, 1, Wv);
    add(x, Wv, 2, x);
    
    //! P = P - W1*W1.transpose();
    transpose(W1, 2, 2, W1t);
    mul(W1, W1t, 2, 2, 2, W1W1t);
    sub(P, W1W1t, 4, P);
}

// Fused Matrix Operations but no AVX
void KF_cholesky_update_fused_ops(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H)
{
    double PHt[4], S[4], St[4], SChol[4], SCholInv[4];
    double W1[4], W[4], W1W1t[4];

    mmT_2x2(P, H, PHt);                 //! PHt = P*Ht
    copy(R, 4, S);                      //! S = R;
    mmadd_2x2(H, PHt, S);               //! S += H*PHt ( S = H*P*H^T + R )

    // optimize or skip
    transpose_2x2(S, St);               //! St = S^T
    add(S, St, 4, S);                   //! S = S + St
    scal(S, 4, 0.5, S);                 //! S = 0.5*S
 
    llt_2x2(S, SChol);                  //! S = SChol * SChol^T
    inv_2x2(SChol, SCholInv);           //! SCholInv = inv( SChol )

#ifdef KF_YGLEE
    mm_2x2(PHt, SCholInv, W1);          //! W1 = PHt * SCholInv
    mmT_2x2(W1, SCholInv, W);
#else
    mmT_2x2(PHt, SCholInv, W1);          //! W1 = PHt * SCholInvt
    mm_2x2(W1, SCholInv, W);
#endif

    mvadd_2x2(W, v, x);                 //! x = x + W*v

    //! P = P - W1*W1.transpose();
    mmT_2x2(W1, W1, W1W1t);
    sub(P, W1W1t, 4, P);
}

// Fused Matrix Operations and AVX
#if __AVX2__
void KF_cholesky_update_fused_ops_avx(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H)
{
    double        S[4] __attribute__ ((aligned(32)));
    double    SChol[4] __attribute__ ((aligned(32)));
    double SCholInv[4] __attribute__ ((aligned(32)));
    
    __m256d p = _mm256_load_pd( P );
    __m256d h = _mm256_load_pd( H );

    __m256d pht = _mmT_2x2_avx_v3(p, h);                        //! PHt = P*Ht
    __m256d s = _mmadd_2x2_avx_v2(h, pht, _mm256_load_pd( R )); //! S += H*PHt ( S = H*P*H^T + R )

    s = _mm256_add_pd( s, _mm256_permute4x64_pd( s, 0b11011000 ) ); // S = S + S^T
    s = _mm256_mul_pd( _mm256_set1_pd( 0.5 ), s );                  // S = 0.5*S

    _mm256_store_pd(S, s);

    // TODO: AVX
    llt_2x2(S, SChol);        //! S = SChol * SChol^T
    inv_2x2(SChol, SCholInv); //! SCholInv = inv( SChol )

    __m256d scholinv = _mm256_load_pd( SCholInv );

#ifdef KF_YGLEE
    __m256d w1 = _mm_2x2_avx_v1(pht, scholinv);  //! W1 = PHt * SCholInv
    __m256d w  = _mmT_2x2_avx_v3(w1, scholinv);
#else
    __m256d w1 = _mmT_2x2_avx_v3(pht, scholinv); //! W1 = PHt * SCholInvt
    __m256d w  = _mm_2x2_avx_v1(w1, scholinv);
#endif

    _mm_store_pd( x, _mvadd_2x2_avx_v1(w, _mm_load_pd(v), _mm_load_pd(x) ) ); //! x = x + W*v

    //! P = P - W1*W1.transpose();
    __m256d w1w1t = _mmT_2x2_avx_v3(w1, w1);
    _mm256_store_pd(P, _mm256_sub_pd( p, w1w1t ));
}
#endif

// ONLY `IF KF_YGLEE` IS NOT DEFINED
// We can substantially reduce the FLOPs by rewritting the function as follows:
/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: 2 neg + 25 adds + 34 muls + 1 divs = 62 flops
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
#ifndef KF_YGLEE
void KF_cholesky_update_reduced_flops(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H)
{
    double PHt[4], S[4], St[4], Sinv[4], W[4], W1W1t[4];

    mmT_2x2(P, H, PHt);
    copy(R, 4, S);        //! S = R;
    mmadd_2x2(H, PHt, S); //! S += H*PHt ( S = H*P*H^T + R )

    // optimize or skip
    //transpose_2x2(S, St); //! St = S^T
    //add(S, St, 4, S);     //! S = S + St
    //scal(S, 4, 0.5, S);   //! S = 0.5*S

    inv_2x2(S, Sinv);     //! Sinv = S^(-1)
    mm_2x2(PHt, Sinv, W); //! W = PHt*Sinv

    mvadd_2x2(W, v, x);   //! x = x + W*v

    mmT_2x2(W, PHt, W1W1t);
    sub(P, W1W1t, 4, P);
}
#endif

#ifndef KF_YGLEE
    #ifdef __AVX2__
void KF_cholesky_update_reduced_flops_avx(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H)
{
    double     S[4] __attribute__ ((aligned(32)));
    double  Sinv[4] __attribute__ ((aligned(32)));

    __m256d p = _mm256_load_pd( P );
    __m256d h = _mm256_load_pd( H );

    __m256d pht = _mmT_2x2_avx_v3(p, h);                        //! PHt = P*Ht
    __m256d s = _mmadd_2x2_avx_v2(h, pht, _mm256_load_pd( R )); //! S += H*PHt ( S = H*P*H^T + R )

    // maybe skip
    //s = _mm256_add_pd( s, _mm256_permute4x64_pd( s, 0b11011000 ) ); // S = S + S^T
    //s = _mm256_mul_pd( _mm256_set1_pd( 0.5 ), s );                  // S = 0.5*S
    _mm256_store_pd(S, s);
 
    inv_2x2(S, Sinv);     //! Sinv = S^(-1) ( TODO: AVX )
    
    __m256d sinv = _mm256_load_pd( Sinv );
    
    __m256d w = _mm_2x2_avx_v1(pht, sinv); //! W = PHt*Sinv

    _mm_store_pd( x, _mvadd_2x2_avx_v1(w, _mm_load_pd(v), _mm_load_pd(x) ) ); //! x = x + W*v

    __m256d w1w1t = _mmT_2x2_avx_v3(w, pht);
    _mm256_store_pd(P, _mm256_sub_pd( p, w1w1t ));
}
    #endif
#endif

// x0x2 holds 2 Vector2d
// x1x3 holds 2 Vector2d
// P0 holds 1 Matrix2d
// P1 holds 1 Matrix2d
// P2 holds 1 Matrix2d
// P3 holds 1 Matrix2d
// v0v2 holds 2 Vector2d
// v1v3 holds 2 Vector2d
// R holds 1 Matrix2d ( same for all consecutive iterations )
// H0 holds 1 Matrix2d
// H1 holds 1 Matrix2d
// H2 holds 1 Matrix2d
// H3 holds 1 Matrix2d
/*
void KF_cholesky_update_unrolled4_avx(__m256d *x0x2,
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

    // TODO: AVX -> Transpose first then compute then back transpose
    inv_2x2(S0, S0inv); //! Sinv = S^(-1)
    inv_2x2(S1, S1inv); //! Sinv = S^(-1)
    inv_2x2(S2, S2inv); //! Sinv = S^(-1)
    inv_2x2(S3, S3inv); //! Sinv = S^(-1)

    //! W = PHt*Sinv
    __m256d const W0 = _mm_2x2_avx_v1(P0H0t, S0inv);
    __m256d const W1 = _mm_2x2_avx_v1(P1H1t, S1inv);
    __m256d const W2 = _mm_2x2_avx_v1(P2H2t, S2inv);
    __m256d const W3 = _mm_2x2_avx_v1(P3H3t, S3inv);

    //! Permutes for stacked mvadds
    __m256d const x2x0 = _mm256_permute2f128_pd( *x0x2, *x0x2, 0b00000001 );
    __m256d const v2v0 = _mm256_permute2f128_pd(  v0v2,  v0v2, 0b00000001 );
    __m256d const x3x1 = _mm256_permute2f128_pd( *x1x3, *x1x3, 0b00000001 );
    __m256d const v3v1 = _mm256_permute2f128_pd(  v1v3,  v1v3, 0b00000001 );

    //! x = x + W*v
    __m128d const x0 =_mvadd_2x2_avx_v1( W0, _mm256_castpd256_pd128(v0v2), _mm256_castpd256_pd128(*x0x2) );
    __m128d const x1 =_mvadd_2x2_avx_v1( W1, _mm256_castpd256_pd128(v1v3), _mm256_castpd256_pd128(*x1x3) );
    __m128d const x2 =_mvadd_2x2_avx_v1( W2, _mm256_castpd256_pd128(v2v0), _mm256_castpd256_pd128( x2x0) );
    __m128d const x3 =_mvadd_2x2_avx_v1( W3, _mm256_castpd256_pd128(v3v1), _mm256_castpd256_pd128( x3x1) );

    //! Store back in x0x2, x1x3
    *x0x2 = _mm256_blend_pd( _mm256_broadcast_pd( &x0 ), _mm256_broadcast_pd( &x2 ), 0b1100 );
    *x1x3 = _mm256_blend_pd( _mm256_broadcast_pd( &x1 ), _mm256_broadcast_pd( &x3 ), 0b1100 );

    //! Z = W * (PHt)^T = P*Ht*Sinv*H*Pt
    __m256d const Z0 = _mmT_2x2_avx_v3( W0, P0H0t );
    __m256d const Z1 = _mmT_2x2_avx_v3( W1, P1H1t );
    __m256d const Z2 = _mmT_2x2_avx_v3( W2, P2H2t );
    __m256d const Z3 = _mmT_2x2_avx_v3( W3, P3H3t );

    _mm256_store_pd( *P0, _mm256_sub_pd( *P0, Z0 ) );
    _mm256_store_pd( *P1, _mm256_sub_pd( *P1, Z1 ) );
    _mm256_store_pd( *P2, _mm256_sub_pd( *P2, Z2 ) );
    _mm256_store_pd( *P3, _mm256_sub_pd( *P3, Z3 ) );
}
*/
