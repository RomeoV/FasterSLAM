#pragma once
#include <iostream>
#include <immintrin.h>

//! ------------------------------------------------------- //
//! ---------- Linear Algebra Utility Functions ----------- //
//! ------------------------------------------------------- //

//HACK: My compiler could not find this function even though its listed in intrinsics.
// Source: https://github.com/gcc-mirror/gcc/blob/master/gcc/config/i386/avxintrin.h
extern __inline void _mm256_store2_m128d (double *__PH, double *__PL, __m256d __A)
{
  _mm_storeu_pd (__PL, _mm256_castpd256_pd128 (__A));
  _mm_storeu_pd (__PH, _mm256_extractf128_pd (__A, 1));
}

//! Prints an array
void print(const double *x, size_t rows, size_t cols, std::ostream& = std::cout);

//! Prints a SIMD vector
void print256d(__m256d var);

//! Fills an array with a specific value
void fill(double *x, size_t size, double val);

//! Copies all values from ref to target
void copy(const double* ref, size_t N, double* target);

//! Fills an array with random values in the range [lo, hi]
void fill_rand(double *x, size_t size, double lo, double hi);

__m256d fill_rand_avx(double lo, double hi);

//! Same as fill_rand but uses a faster pseudo RNG
void fill_rand_fast(double *x, size_t size, double lo, double hi);

//! ------------------------------------------------------- //
//! -------------- Basic Matrix Operations ---------------- //
//! ------------------------------------------------------- //

//! Matrix Transpose
//! Matrix A ( mA x nA ) -> Input
//! Matrix T ( nA x mA ) -> Output
//! CAUTION!!! For in place transposition use stranspose() instead
void transpose(const double *A, size_t mA, size_t nA, double *T); 

//! Matrix Transpose ( 2x2 )
void transpose_2x2(const double *A, double *T);

//! Matrix Self Transpose ( 2x2, in place )
void stranspose_2x2(double *A);

__m256d _transpose_2x2_avx( __m256d A );

//! Adds two arrays
void add(const double *x, const double *y, size_t size, double* z); 

//! Subtracts two arrays
void sub(const double *x, const double *y, size_t size, double* z);

//! Scales an array by a scalar
void scal(const double *x, size_t size, double a, double *y);

//! Matrix x Matrix Multiplication: C = A * B 
//! Matrix A ( mA x nA ) -> Input
//! Matrix B ( nA x nB ) -> Input
//! Matrix C ( mA x nB ) -> Output
void mul(const double *A, const double *B, size_t mA, size_t nA, size_t nB, double *C);

//! Matrix x Matrix Multiplication ( 2x2 )
void mm_2x2(const double *A, const double *B, double *C);
void mm_2x2_avx_v1(const double *A, const double *B, double *C);
void mm_2x2_avx_v2(const double *A, const double *B, double *C);
void mm_2x2_avx_v3(const double *A, const double *B, double *C);
__m256d _mm_2x2_avx_v1( __m256d a, __m256d b );

//! Matrix x Matrix Transpose Multiplication ( 2x2 )
void mmT_2x2(const double *A, const double *B, double *C);
void mmT_2x2_avx_v1(const double *A, const double *B, double *C);
void mmT_2x2_avx_v2(const double *A, const double *B, double *C);
void mmT_2x2_avx_v3(const double *A, const double *B, double *C);
__m256d _mmT_2x2_avx_v3( __m256d a0123, __m256d b0123 );

//! C += A*B ( 2x2 )
void mmadd_2x2(const double *A, const double *B, double *C);
void mmadd_2x2_avx_v1(const double *A, const double *B, double *C);
void mmadd_2x2_avx_v2(const double *A, const double *B, double *C);
__m256d _mmadd_2x2_avx_v2(__m256d a, __m256d b, __m256d c);

__m256d _mmTadd_2x2_avx_v2(__m256d a0123, __m256d b0123, __m256d c);

//! Matrix x Vector Multiplication ( 2x2 )
void mv_2x2(const double *A, const double *b, double *c);
__m128d _mv_2x2_avx_v1( __m256d const a, __m128d const b );

//! c += A*b ( 2x2 )
void mvadd_2x2(const double *A, const double *b, double *c);
__m128d _mvadd_2x2_avx_v1( __m256d const a, __m128d const b, __m128d const c );

//! Cholesky Factorization of a 2x2 SPD Matrix A = L * L^T, L lower triangular
void llt_2x2(const double *A, double *L);

//! Inverse of a 2x2 Matrix
void inv_2x2(const double *A, double *Ainv);

//! Determinant for 2x2 matrix
//! @param A pointer to row major continuous memory
double determinant_2x2(const double* A);


// Instrumenting flops and memory count

double copy_flops(const double* ref, size_t N, double* target);
double copy_memory(const double* ref, size_t N, double* target);

double fill_rand_flops(double *x, size_t size, double lo, double hi);
double fill_rand_memory(double *x, size_t size, double lo, double hi);

double transpose_flops(const double *A, size_t mA, size_t nA, double *T);
double transpose_memory(const double *A, size_t mA, size_t nA, double *T);


double add_flops(const double *x, const double *y, size_t size, double* z);
double add_memory(const double *x, const double *y, size_t size, double* z);


double mul_flops(const double *A, const double *B, size_t mA, size_t nA, size_t nB, double *C);
double mul_memory(const double *A, const double *B, size_t mA, size_t nA, size_t nB, double *C);

double llt_2x2_flops(const double *A, double *L);
double llt_2x2_memory(const double *A, double *L);

double inv_2x2_flops(const double *A, double *Ainv);
double inv_2x2_memory(const double *A, double *Ainv);

double determinant_2x2_flops(const double* A);
double determinant_2x2_memory(const double* A);

/** Executes four v.T @ M @ v.T calculations
 * @param m1 - m4 Row major matrices
 * @param v12 vector v1 and v2 in order
 * @param v34 vector 3 and v4 in order
 *
 * @result AVX vector with the results
 */
__m256d mm_vT_M_v_avx2(const __m256d m1,  const __m256d m2,
                       const __m256d m3,  const __m256d m4,
                       const __m256d v12, const __m256d v34);

__m256d mm_vT_M_v_avx2_phil(const __m256d m1,  const __m256d m2,
                       const __m256d m3,  const __m256d m4,
                       const __m256d v12, const __m256d v34);
//! Takes four registers and returns their transpose equivalents
void register_transpose(__m256d const r0,
                        __m256d const r1,
                        __m256d const r2,
                        __m256d const r3,
                        __m256d *t0, __m256d *t1, __m256d *t2, __m256d *t3);

//! Takes four registers each of which holds a matrix 2x2 (row-major),
//! and returns their inverses
void batch_inverse_2x2(__m256d const r0,
                       __m256d const r1,
                       __m256d const r2,
                       __m256d const r3,
                       __m256d *inv0,
                       __m256d *inv1,
                       __m256d *inv2,
                       __m256d *inv3);

