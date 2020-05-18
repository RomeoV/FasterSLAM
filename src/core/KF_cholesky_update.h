#pragma once
#include "typedefs.h"
#include <immintrin.h>

/*!
    Update Kalman Filter using Cholesky decomposition. 
    @param[out] x 	    State vector.
    @param[out] P 	    Covariance matrix of state.
    @param[in]  v   	Velocity at current state.
    @param[in]  R    	Covariance matrix of measurements.
    @param[in]  H    	Jacobian of h wrt feature states.
*/
void KF_cholesky_update(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);

void KF_cholesky_update_base(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);
void KF_cholesky_update_fused_ops(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);
void KF_cholesky_update_fused_ops_avx(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);

#ifndef KF_YGLEE
void KF_cholesky_update_reduced_flops(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);
void KF_cholesky_update_reduced_flops_avx(Vector2d x, Matrix2d P, cVector2d v, cMatrix2d R, cMatrix2d H);

void KF_cholesky_update_unrolled4_avx(__m256d *x0x2,
                                      __m256d *x1x3,
                                      __m256d *P0,
                                      __m256d *P1,
                                      __m256d *P2,
                                      __m256d *P3,
                                      __m256d const v0v2,
                                      __m256d const v1v3,
                                      __m256d const R,
                                      __m256d const H0,
                                      __m256d const H1,
                                      __m256d const H2,
                                      __m256d const H3);
#endif
