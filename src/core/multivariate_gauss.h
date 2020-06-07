#pragma once
#include "typedefs.h"

/*!
    Computes Multivariate Gaussian for ONE sample.
    @param[in]  x        Mean vector.
    @param[in]  P        Covariance matrix.
    @param[out] result   Sample set.
 */

void multivariate_gauss(cVector2d x, cMatrix2d P, Vector2d result);

void multivariate_gauss_base(cVector2d x, cMatrix2d P, Vector2d result);

void multivariate_gauss_active(cVector2d x, cMatrix2d P, Vector2d result);

FlopCount multivariate_gauss_base_flops(cVector2d x, cMatrix2d P, Vector2d result);
double multivariate_gauss_base_memory(cVector2d x, cMatrix2d P, Vector2d result);

FlopCount multivariate_gauss_active_flops(cVector2d x, cMatrix2d P, Vector2d result);
double multivariate_gauss_active_memory(cVector2d x, cMatrix2d P, Vector2d result);


