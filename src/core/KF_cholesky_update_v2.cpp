#include "KF_cholesky_update.h"
#include "linalg.h"
#include <immintrin.h>

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

void KF_cholesky_update(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H) {
    //KF_cholesky_update_base(x, P, v, R, H);
    //KF_cholesky_update_v1(x, P, v, R, H);
#ifdef __AVX__
    //KF_cholesky_update_v2_avx_v1(x, P, v, R, H);
    KF_cholesky_update_v2_avx_v2(x, P, v, R, H);
#else
    KF_cholesky_update_v2(x, P, v, R, H);
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

    transpose(SCholInv, 2, 2, SCholInvt); //! SCholInvt = SCholInv.transpose();
    mul(PHt, SCholInvt, 2, 2, 2, W1);     //! W1 = PHt * SCholInvt;
    mul(W1, SCholInv, 2, 2, 2, W);        //! W = W1 * SCholInv;

    //! x = x + Wv;
    mul(W, v, 2, 2, 1, Wv);
    add(x, Wv, 2, x);
    
    //! P = P - W1*W1.transpose();
    transpose(W1, 2, 2, W1t);
    mul(W1, W1t, 2, 2, 2, W1W1t);
    sub(P, W1W1t, 4, P);
}

// maybe more stable (according to MATLAB code)
// TODO count flops
void KF_cholesky_update_v1(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H)
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

    mmT_2x2(PHt, SCholInv, W1);          //! W1 = PHt * SCholInv
    mm_2x2(W1, SCholInv, W);

    mvadd_2x2(W, v, x);                 //! x = x + W*v

    //! P = P - W1*W1.transpose();
    mmT_2x2(W1, W1, W1W1t);
    sub(P, W1W1t, 4, P);
}

// fast
/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: 2 neg + 25 adds + 34 muls + 1 divs = 62 flops
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void KF_cholesky_update_v2(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H)
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

#ifdef __AVX__
void KF_cholesky_update_v2_avx_v1(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H)
{
    double   PHt[4] __attribute__ ((aligned(32)));
    double     S[4] __attribute__ ((aligned(32)));
    double    St[4] __attribute__ ((aligned(32)));
    double  Sinv[4] __attribute__ ((aligned(32)));
    double     W[4] __attribute__ ((aligned(32)));
    double W1W1t[4] __attribute__ ((aligned(32)));

    // ------------------- //
    // mmT_2x2(P, H, PHt); //
    // ------------------- //
    __m256d p = _mm256_load_pd( P );
    __m256d h = _mm256_load_pd( H );

    __m256d pht = _mm256_mul_pd(
                  _mm256_permute_pd( p, 0b0000 ),
                  _mm256_permute4x64_pd( h, 0b10001000 ) );

#ifdef __FMA__
    pht = _mm256_fmadd_pd(
           _mm256_permute_pd( p, 0b1111 ),
           _mm256_permute4x64_pd( h, 0b11011101 ), pht );
#else
    pht = _mm256_add_pd(
            _mm256_mul_pd(
               _mm256_permute_pd( p, 0b1111 ),
               _mm256_permute4x64_pd( h, 0b11011101 )
           ), pht );
#endif
 
    // --------------------- //
    // copy(R, 4, S);        //! S = R;
    // mmadd_2x2(H, PHt, S); //! S += H*PHt ( S = H*P*H^T + R )
    // --------------------- //
    __m256d s = _mm256_load_pd( R );
    
#ifdef __FMA__
    s = _mm256_fmadd_pd(
           _mm256_permute_pd( h, 0b0000 ),
           _mm256_permute2f128_pd( pht, pht, 0b00000000 ), s );
#else
    s = _mm256_add_pd(
            _mm256_mul_pd(
               _mm256_permute_pd( h, 0b0000 ),
               _mm256_permute2f128_pd( pht, pht, 0b00000000 )
           ), s );
#endif

#ifdef __FMA__
    s = _mm256_fmadd_pd(
           _mm256_permute_pd( h, 0b1111 ),
           _mm256_permute2f128_pd( pht, pht, 0b01010101 ), s );
#else
    s = _mm256_add_pd(
            _mm256_mul_pd(
               _mm256_permute_pd( h, 0b1111 ),
               _mm256_permute2f128_pd( pht, pht, 0b01010101 )
           ), s );
#endif

    _mm256_store_pd( S, s );

    // optimize or skip
    //transpose_2x2(S, St); //! St = S^T
    //add(S, St, 4, S);     //! S = S + St
    //scal(S, 4, 0.5, S);   //! S = 0.5*S

    inv_2x2(S, Sinv);     //! Sinv = S^(-1) ( TODO: AVX )

    // --------------------- //
    // mm_2x2(PHt, Sinv, W); //! W = PHt*Sinv
    // --------------------- //
    
    __m256d sinv = _mm256_load_pd( Sinv );

    __m256d w = _mm256_mul_pd(
                   _mm256_permute_pd( pht, 0b0000 ),
                   _mm256_permute2f128_pd( sinv, sinv, 0b00000000 ) );

#ifdef __FMA__
    w = _mm256_fmadd_pd(
           _mm256_permute_pd( pht, 0b1111 ),
           _mm256_permute2f128_pd( sinv, sinv, 0b01010101 ), w );
#else
    w = _mm256_add_pd(
            _mm256_mul_pd(
               _mm256_permute_pd( pht, 0b1111 ),
               _mm256_permute2f128_pd( sinv, sinv, 0b01010101 )
           ), w );
#endif

    _mm256_store_pd( W, w );    

    mvadd_2x2(W, v, x);   //! x = x + W*v ( TODO: AVX )
 
    // ----------------------- //
    // mmT_2x2(W, PHt, W1W1t); //
    // ----------------------- //
    __m256d w1w1t = _mm256_mul_pd(
                       _mm256_permute_pd( w, 0b0000 ),
                       _mm256_permute4x64_pd( pht, 0b10001000 ) );

#ifdef __FMA__
    w1w1t = _mm256_fmadd_pd(
               _mm256_permute_pd( w, 0b1111 ),
               _mm256_permute4x64_pd( pht, 0b11011101 ), w1w1t );
#else
    w1w1t = _mm256_add_pd(
            _mm256_mul_pd(
                   _mm256_permute_pd( w, 0b1111 ),
                   _mm256_permute4x64_pd( pht, 0b11011101 )
               ), w1w1t );
#endif
    
    // -------------------- //
    // sub(P, W1W1t, 4, P); //
    // -------------------- //
    _mm256_store_pd(P, _mm256_sub_pd( p, w1w1t ));
}
#endif

#ifdef __AVX__
// fast
void KF_cholesky_update_v2_avx_v2(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H)
{
    double   PHt[4] __attribute__ ((aligned(32)));
    double     S[4] __attribute__ ((aligned(32)));
    double    St[4] __attribute__ ((aligned(32)));
    double  Sinv[4] __attribute__ ((aligned(32)));
    double     W[4] __attribute__ ((aligned(32)));
    double W1W1t[4] __attribute__ ((aligned(32)));

    // mmT_2x2(P, H, PHt);
    __m256d p = _mm256_load_pd( P ); // should be aligned
    __m256d h = _mm256_load_pd( H ); // should be aligned
    __m256d p0022 = _mm256_permute_pd( p, 0b0000 );
    __m256d p1133 = _mm256_permute_pd( p, 0b1111 );
    __m256d h0202 = _mm256_permute4x64_pd( h, 0b10001000 );
    __m256d h1313 = _mm256_permute4x64_pd( h, 0b11011101 );
#ifdef __FMA__
    __m256d pht = _mm256_fmadd_pd( p1133, h1313, _mm256_mul_pd( p0022, h0202 ) );
#else
    __m256d pht = _mm256_add_pd( _m256_mul_pd(p1133, h1313), _mm256_mul_pd( p0022, h0202 ) );
#endif

    // copy(R, 4, S);        //! S = R; 
    __m256d s = _mm256_load_pd( R ); // should be aligned
    
    // mmadd_2x2(H, PHt, S); //! S += H*PHt ( S = H*P*H^T + R )
    __m256d h0022 = _mm256_permute_pd( h, 0b0000 );
    __m256d h1133 = _mm256_permute_pd( h, 0b1111 );
    __m256d pht0101 = _mm256_permute2f128_pd( pht, pht, 0b00000000 );
    __m256d pht2323 = _mm256_permute2f128_pd( pht, pht, 0b01010101 );
#ifdef __FMA__
    s = _mm256_fmadd_pd( h0022, pht0101, s );
#else
    s = _mm256_add_pd( _m256_mul_pd(h0022, pht0101), s );
#endif
#ifdef __FMA__
    s = _mm256_fmadd_pd( h1133, pht2323, s );
#else
    s = _mm256_add_pd( _m256_mul_pd(h1133, pht2323), s );
#endif
    _mm256_store_pd( S, s );

    inv_2x2(S, Sinv);     //! Sinv = S^(-1) ( TODO: AVX )
     
    // mm_2x2(PHt, Sinv, W); //! W = PHt*Sinv
    __m256d sinv = _mm256_load_pd( Sinv );
    __m256d pht0022 = _mm256_permute_pd( pht, 0b0000 );
    __m256d pht1133 = _mm256_permute_pd( pht, 0b1111 );
    __m256d sinv0101 = _mm256_permute2f128_pd( sinv, sinv, 0b00000000 );
    __m256d sinv2323 = _mm256_permute2f128_pd( sinv, sinv, 0b01010101 );
#ifdef __FMA__
    __m256d w = _mm256_fmadd_pd( pht1133, sinv2323, _mm256_mul_pd( pht0022, sinv0101 ) ); 
#else
    __m256d w = _mm256_add_pd( _m256_mul_pd(pht1133, sinv2323), _mm256_mul_pd( pht0022, sinv0101 ) ); 
#endif
    _mm256_store_pd(W, w);

    mvadd_2x2(W, v, x);   //! x = x + W*v ( TODO: AVX )

    // mmT_2x2(W, PHt, W1W1t); 
    __m256d w0022 = _mm256_permute_pd( w, 0b0000 );
    __m256d w1133 = _mm256_permute_pd( w, 0b1111 );
    __m256d pht0202 = _mm256_permute4x64_pd( pht, 0b10001000 );
    __m256d pht1313 = _mm256_permute4x64_pd( pht, 0b11011101 ); 
#ifdef __FMA__
    __m256d w1w1t = _mm256_fmadd_pd( w1133, pht1313, _mm256_mul_pd( w0022, pht0202 ) );
#else
    __m256d w1w1t = _mm256_add_pd( _m256_mul_pd(w1133, pht1313), _mm256_mul_pd( w0022, pht0202 ) );
#endif
 
    // sub(P, W1W1t, 4, P);
    _mm256_store_pd(P, _mm256_sub_pd( p, w1w1t ));
}
#endif
