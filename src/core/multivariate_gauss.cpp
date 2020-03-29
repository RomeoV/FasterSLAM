#include "multivariate_gauss.h"
#include "linalg.h"

void multivariate_gauss(Vector2d x, Matrix2d P, Vector2d result)
{
    double S[4]; // S is a 2x2 Matrix
    llt_2x2(P, S); // P = S x S^T, S lower triangular

    double X[2];  
    fill_rand(X, 2, -1.0, 1.0);
    
    // result = S*X + x
    mul(S, X, 2, 2, 1, result);
    add(result, x, 2, result);
}
