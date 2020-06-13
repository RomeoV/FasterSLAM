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

void KF_cholesky_update_active(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H) {
    KF_cholesky_update(x, P, v, R, H);
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


// Utils

// Naive base flops
double KF_cholesky_update_base_flops(Vector2d x, Matrix2d P,
        cVector2d v, cMatrix2d R, cMatrix2d H) {
    double Ht[4], PHt[4], HPHt[4];
    double S[4], St[4], SChol[4], SCholInv[4], SCholInvt[4];
    double W1[4], W[4], Wv[2], W1t[4], W1W1t[4];
    double flops = 0.0;

    flops += transpose_flops(H, 2, 2, Ht);               //! Ht = H.transpose();
    flops += mul_flops(P, Ht, 2, 2, 2, PHt);             //! PHt = P*Ht;
    flops += mul_flops(H, PHt, 2, 2, 2, HPHt);           //! HPHt = H*PHt;
    flops += add_flops(HPHt, R, 4, S);                   //! S = HPHt + R;

    flops += transpose_flops(S, 2, 2, St);               //! St = S.transpose();
    flops += add_flops(S, St, 4, S);                     //! S = S + St
    flops += scal_flops(S, 4, 0.5, S);                   //! S = 0.5*S;

    flops += llt_2x2_flops(S, SChol);                    //! S = SChol * SChol^T
    flops += inv_2x2_flops(SChol, SCholInv);             //! SCholInv = inv( SChol )

#ifdef KF_YGLEE
    flops += mul_flops(PHt, SCholInv, 2, 2, 2, W1);     //! W1 = PHt * SCholInv;
    flops += transpose_flops(SCholInv, 2, 2, SCholInvt); //! SCholInvt = SCholInv.transpose();
    flops += mul_flops(W1, SCholInvt, 2, 2, 2, W);        //! W = W1 * SCholInvt;
#else
    flops += transpose_flops(SCholInv, 2, 2, SCholInvt); //! SCholInvt = SCholInv.transpose();
    flops += mul_flops(PHt, SCholInvt, 2, 2, 2, W1);     //! W1 = PHt * SCholInvt;
    flops += mul_flops(W1, SCholInv, 2, 2, 2, W);        //! W = W1 * SCholInv;
#endif

    //! x = x + Wv;
    flops += mul_flops(W, v, 2, 2, 1, Wv);
    flops += add_flops(x, Wv, 2, x);
 
    //! P = P - W1*W1.transpose();
    flops += transpose_flops(W1, 2, 2, W1t);
    flops += mul_flops(W1, W1t, 2, 2, 2, W1W1t);
    flops += sub_flops(P, W1W1t, 4, P);

    return flops;
}

// Naive base memory
double KF_cholesky_update_base_memory(Vector2d x, Matrix2d P,
        cVector2d v, cMatrix2d R, cMatrix2d H) {
    double Ht[4], PHt[4], HPHt[4];
    double S[4], St[4], SChol[4], SCholInv[4], SCholInvt[4];
    double W1[4], W[4], Wv[2], W1t[4], W1W1t[4];
    double memory = 0.0;

    memory += transpose_memory(H, 2, 2, Ht);               //! Ht = H.transpose();
    memory += mul_memory(P, Ht, 2, 2, 2, PHt);             //! PHt = P*Ht;
    memory += mul_memory(H, PHt, 2, 2, 2, HPHt);           //! HPHt = H*PHt;
    memory += add_memory(HPHt, R, 4, S);                   //! S = HPHt + R;

    memory += transpose_memory(S, 2, 2, St);               //! St = S.transpose();
    memory += add_memory(S, St, 4, S);                     //! S = S + St
    memory += scal_memory(S, 4, 0.5, S);                   //! S = 0.5*S;

    memory += llt_2x2_memory(S, SChol);                    //! S = SChol * SChol^T
    memory += inv_2x2_memory(SChol, SCholInv);             //! SCholInv = inv( SChol )

#ifdef KF_YGLEE
    memory += mul_memory(PHt, SCholInv, 2, 2, 2, W1);     //! W1 = PHt * SCholInv;
    memory += transpose_memory(SCholInv, 2, 2, SCholInvt); //! SCholInvt = SCholInv.transpose();
    memory += mul_memory(W1, SCholInvt, 2, 2, 2, W);        //! W = W1 * SCholInvt;
#else
    memory += transpose_memory(SCholInv, 2, 2, SCholInvt); //! SCholInvt = SCholInv.transpose();
    memory += mul_memory(PHt, SCholInvt, 2, 2, 2, W1);     //! W1 = PHt * SCholInvt;
    memory += mul_memory(W1, SCholInv, 2, 2, 2, W);        //! W = W1 * SCholInv;
#endif

    //! x = x + Wv;
    memory += mul_memory(W, v, 2, 2, 1, Wv);
    memory += add_memory(x, Wv, 2, x);
 
    //! P = P - W1*W1.transpose();
    memory += transpose_memory(W1, 2, 2, W1t);
    memory += mul_memory(W1, W1t, 2, 2, 2, W1W1t);
    memory += sub_memory(P, W1W1t, 4, P);

    return 3*2+3*4; // only inputs/outputs ( one read + write each )
    return memory;
}

// Active flops
double KF_cholesky_update_active_flops(Vector2d x, Matrix2d P,
        cVector2d v, cMatrix2d R, cMatrix2d H) {

#ifdef KF_YGLEE
    // SHOULD !!NOT!! BE USED
    // Just fusing some ops
    double PHt[4], S[4], St[4], SChol[4], SCholInv[4];
    double W1[4], W[4], W1W1t[4];
    double flops = 0.0;

    flops += mmT_2x2_flops(P, H, PHt);                 //! PHt = P*Ht
    flops += copy_flops(R, 4, S);                      //! S = R;
    flops += mmadd_2x2_flops(H, PHt, S);               //! S += H*PHt ( S = H*P*H^T + R )

    // optimize or skip
    flops += transpose_2x2_flops(S, St);               //! St = S^T
    flops += add_flops(S, St, 4, S);                   //! S = S + St
    flops += scal_flops(S, 4, 0.5, S);                 //! S = 0.5*S
 
    flops += llt_2x2_flops(S, SChol);                  //! S = SChol * SChol^T
    flops += inv_2x2_flops(SChol, SCholInv);           //! SCholInv = inv( SChol )

    flops += mm_2x2_flops(PHt, SCholInv, W1);          //! W1 = PHt * SCholInv
    flops += mmT_2x2_flops(W1, SCholInv, W);

    flops += mvadd_2x2_flops(W, v, x);                 //! x = x + W*v

    //! P = P - W1*W1.transpose();
    flops += mmT_2x2_flops(W1, W1, W1W1t);
    flops += sub_flops(P, W1W1t, 4, P);

    return flops;

#else    
    // THIS IS WHAT SHOULD BE USED
    // After some math reformulation and fusing some linalg ops
    // This is based on KF_cholesky_update_reduced_flops
    double PHt[4], S[4], St[4], Sinv[4], W[4], W1W1t[4];
    double flops = 0.0;

    flops += mmT_2x2_flops(P, H, PHt);
    flops += copy_flops(R, 4, S);        //! S = R;
    flops += mmadd_2x2_flops(H, PHt, S); //! S += H*PHt ( S = H*P*H^T + R )

    // optimize or skip
    //flops += transpose_2x2_flops(S, St); //! St = S^T
    //flops += add_flops(S, St, 4, S);     //! S = S + St
    //flops += scal_flops(S, 4, 0.5, S);   //! S = 0.5*S

    flops += inv_2x2_flops(S, Sinv);     //! Sinv = S^(-1)
    flops += mm_2x2_flops(PHt, Sinv, W); //! W = PHt*Sinv

    flops += mvadd_2x2_flops(W, v, x);   //! x = x + W*v

    flops += mmT_2x2_flops(W, PHt, W1W1t);
    flops += sub_flops(P, W1W1t, 4, P);

    return flops;

#endif
}

// Active memory
double KF_cholesky_update_active_memory(Vector2d x, Matrix2d P,
        cVector2d v, cMatrix2d R, cMatrix2d H) {

#ifdef KF_YGLEE
    // SHOULD !!NOT!! BE USED
    // Just fusing some ops
    double PHt[4], S[4], St[4], SChol[4], SCholInv[4];
    double W1[4], W[4], W1W1t[4];
    double memory = 0.0;

    memory += mmT_2x2_memory(P, H, PHt);                 //! PHt = P*Ht
    memory += copy_memory(R, 4, S);                      //! S = R;
    memory += mmadd_2x2_memory(H, PHt, S);               //! S += H*PHt ( S = H*P*H^T + R )

    // optimize or skip
    memory += transpose_2x2_memory(S, St);               //! St = S^T
    memory += add_memory(S, St, 4, S);                   //! S = S + St
    memory += scal_memory(S, 4, 0.5, S);                 //! S = 0.5*S
 
    memory += llt_2x2_memory(S, SChol);                  //! S = SChol * SChol^T
    memory += inv_2x2_memory(SChol, SCholInv);           //! SCholInv = inv( SChol )

    memory += mm_2x2_memory(PHt, SCholInv, W1);          //! W1 = PHt * SCholInv
    memory += mmT_2x2_memory(W1, SCholInv, W);

    memory += mvadd_2x2_memory(W, v, x);                 //! x = x + W*v

    //! P = P - W1*W1.transpose();
    memory += mmT_2x2_memory(W1, W1, W1W1t);
    memory += sub_memory(P, W1W1t, 4, P);

    return memory;

#else    
    // THIS IS WHAT SHOULD BE USED
    // After some math reformulation and fusing some linalg ops
    // This is based on KF_cholesky_update_reduced_memory
    double PHt[4], S[4], St[4], Sinv[4], W[4], W1W1t[4];
    double memory = 0.0;

    memory += mmT_2x2_memory(P, H, PHt);
    memory += copy_memory(R, 4, S);        //! S = R;
    memory += mmadd_2x2_memory(H, PHt, S); //! S += H*PHt ( S = H*P*H^T + R )

    // optimize or skip
    //memory += transpose_2x2_memory(S, St); //! St = S^T
    //memory += add_memory(S, St, 4, S);     //! S = S + St
    //memory += scal_memory(S, 4, 0.5, S);   //! S = 0.5*S

    memory += inv_2x2_memory(S, Sinv);     //! Sinv = S^(-1)
    memory += mm_2x2_memory(PHt, Sinv, W); //! W = PHt*Sinv

    memory += mvadd_2x2_memory(W, v, x);   //! x = x + W*v

    memory += mmT_2x2_memory(W, PHt, W1W1t);
    memory += sub_memory(P, W1W1t, 4, P);

    return memory;

#endif
}
