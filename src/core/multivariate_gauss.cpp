#include "multivariate_gauss.h"
#include "linalg.h"
#include <math.h>
/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/


void multivariate_gauss(cVector2d x, cMatrix2d P, Vector2d result) {
    multivariate_gauss_active(x, P, result);
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


void multivariate_gauss_fast(cVector2d x, cMatrix2d P, Vector2d result)
{
    double S[4]; //! 2x2 matrix, lower triangular cholesky factor
    //llt_2x2(P, S); //! P = S * S^T

    S[0] = sqrt( P[0] );             // 0*2+0 -> (0,0)
    S[1] = 0.0;                      // 0*2+1 -> (0,1)
    S[2] = P[2] / S[0];              // 1*2+0 -> (1,0)
    S[3] = sqrt( P[3] - S[2]*S[2] ); // 1*2+1 -> (1,1)

    double X[2]; //! 2-vector
    fill_rand(X, 2, -1.0, 1.0);
    
    //! result = S*X + x
    mul(S, X, 2, 2, 1, result);
    add(result, x, 2, result);
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
void multivariate_gauss_active(cVector2d x, cMatrix2d P, Vector2d result)
{
    double S[4]; //! 2x2 matrix, lower triangular cholesky factor
    llt_2x2(P, S); //! P = S * S^T

    double X[2]; //! 2-vector
    fill_rand(X, 2, -1.0, 1.0);
    
    //! result = S*X + x
    result[0] = x[0];
    result[1] = x[1];
    mvadd_2x2(S, X, result);
}
