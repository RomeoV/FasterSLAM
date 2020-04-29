#include "multivariate_gauss.h"
#include "linalg.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/


void multivariate_gauss(cVector2d x, cMatrix2d P, Vector2d result) {
    multivariate_gauss_base(x, P, result);
}

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: 2 sqrt + 2 div + 7 add + 5 mult + 1 neg  = 17 flops
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void multivariate_gauss_base(cVector2d x, cMatrix2d P, Vector2d result)
{
    double S[4]; //! 2x2 matrix, lower triangular cholesky factor
    llt_2x2(P, S); //! P = S * S^T

    double X[2]; //! 2-vector
    fill_rand(X, 2, -1.0, 1.0);
    
    //! result = S*X + x
    mul(S, X, 2, 2, 1, result);
    add(result, x, 2, result);
}
