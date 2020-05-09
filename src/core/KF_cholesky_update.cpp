#include "KF_cholesky_update.h"
#include "linalg.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: 2 neg + 60 adds + 54 muls + 2 divs + 2 sqrts = 120 flops
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/

void KF_cholesky_update(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H) {
    //KF_cholesky_update_base(x, P, v, R, H);
    KF_cholesky_update_v1(x, P, v, R, H);
    //KF_cholesky_update_v2(x, P, v, R, H);
}

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

    mul(PHt, SCholInv, 2, 2, 2, W1);      //! W1 = PHt * SCholInv;
    transpose(SCholInv, 2, 2, SCholInvt); //! SCholInvt = SCholInv.transpose();
    mul(W1, SCholInvt, 2, 2, 2, W);       //! W = W1 * SCholInvt;

    //! x = x + Wv;
    mul(W, v, 2, 2, 1, Wv);
    add(x, Wv, 2, x);
    
    //! P = P - W1*W1.transpose();
    transpose(W1, 2, 2, W1t);
    mul(W1, W1t, 2, 2, 2, W1W1t);
    sub(P, W1W1t, 4, P);
}

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

    mm_2x2(PHt, SCholInv, W1);          //! W1 = PHt * SCholInv
    mmT_2x2(W1, SCholInv, W);

    mvadd_2x2(W, v, x);                 //! x = x + W*v

    //! P = P - W1*W1.transpose();
    mmT_2x2(W1, W1, W1W1t);
    sub(P, W1W1t, 4, P);
}

// slightly modified version ( fails test )
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
