#include "KF_cholesky_update.h"
#include "linalg.h"

void KF_cholesky_update(Vector2d x, Matrix2d P, Vector2d v, Matrix2d R, Matrix2d H)
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
    scal(S, 4, 0.5);                      //! S = 0.5*S;

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
