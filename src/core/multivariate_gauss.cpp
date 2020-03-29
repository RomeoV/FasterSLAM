#include "multivariate_gauss.h"
#include "linalg.h"

void multivariate_gauss(Vector2d x, Matrix2d P, Vector2d result)
{
    double S[4]; //! 2x2 matrix, lower triangular cholesky factor
    llt_2x2(P, S); //! P = S * S^T

    double X[2]; //! 2-vector
    fill_rand(X, 2, -1.0, 1.0);
    
    //! result = S*X + x
    mul(S, X, 2, 2, 1, result);
    add(result, x, 2, result);
}
