#include "KF_cholesky_update.h"

void KF_cholesky_update(Vector2d &x, Matrix2d &P, Vector2d v, Matrix2d R, Matrix2d H)
{
    Matrix2d PHt = P*H.transpose();
    Matrix2d S = H*PHt + R;
    
    S = (S+S.transpose()) * 0.5; //make symmetric
    Matrix2d SChol = S.llt().matrixL();
    SChol.transpose();
    SChol.conjugate();

    Matrix2d SCholInv = SChol.inverse(); //tri matrix
    Matrix2d W1 = PHt * SCholInv;
    Matrix2d W = W1 * SCholInv.transpose();

    x = x + W*v;
    P = P - W1*W1.transpose();
}
