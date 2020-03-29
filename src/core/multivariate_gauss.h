#pragma once
#include "typedefs.h"

/*!
    Computes Multivariate Gaussian for ONE sample.
    @param[out] result   Sample set.
    @param[in]  x        Mean vector.
    @param[in]  P        Covariance matrix.
 */

void multivariate_gauss(Vector2d x, Matrix2d P, Vector2d result);
