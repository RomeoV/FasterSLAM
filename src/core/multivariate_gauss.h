#pragma once
#include "typedefs.h"

/*!
    Computes Multivariate Gaussian for ONE sample.
    @param[out] result   Sample set.
    @param[in]  x        Mean vector.
    @param[in]  P        Covariance matrix.
 */

void multivariate_gauss(cVector2d x, cMatrix2d P, Vector2d result);

void multivariate_gauss_base(cVector2d x, cMatrix2d P, Vector2d result);

void multivariate_gauss_fast_rand(cVector2d x, cMatrix2d P, Vector2d result);