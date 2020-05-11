#include "linalg.h"
#include <iostream>
#include <immintrin.h>
#include <cmath>
#include <cassert>

//! ------------------------------------------------------- //
//! ---------- Linear Algebra Utility Functions ----------- //
//! ------------------------------------------------------- //


/////////////////////////////////////////////////////////////
//! I did not change these functions to base since this would
//! overcomplify the project. Do not edit any functions here, 
//! rather add new ones if necessary, or do it in place!!!
/////////////////////////////////////////////////////////////


//! Prints an array
void print(const double *x, size_t rows, size_t cols, std::ostream& stream) {
    assert( x != NULL );
    for (size_t i = 0; i < rows; i++) {
        stream << std::endl;
        for (size_t j = 0; j < cols; j++) {
            stream << "\t" << x[i*cols+j];
        }
    }
    stream << std::endl;
}

//! Fills an array with a specific value
void fill(double *x, size_t size, double val) {
    assert( x != NULL );
    for (size_t i = 0; i < size; i++) {
        x[i] = val; 
    }
}

//! Copies all values from ref to target
void copy(const double* ref, size_t N, double* target) {
    for (size_t i = 0; i < N; i++) {
        target[i] = ref[i];
    }
}

//! Fills an array with random values in the range [lo, hi]
void fill_rand(double *x, size_t size, double lo, double hi) {
    assert( x != NULL && hi > lo );
    double range = hi - lo;
    for (size_t i = 0; i < size; i++) {
        x[i] = lo + range * ((double)rand())/((double)RAND_MAX);
    }
}

//! ------------------------------------------------------- //
//! -------------- Basic Matrix Operations ---------------- //
//! ------------------------------------------------------- //

//! Matrix Transpose
void transpose(const double *A, size_t mA, size_t nA, double *T) {
    assert( A != NULL && T != NULL );
    for (size_t i = 0; i < nA; i++) {
        for (size_t j = 0; j < mA; j++) {
             T[i*mA+j] = A[j*nA+i];
        }
    }
}

void transpose_2x2(const double *A, double *T) {
    T[0] = A[0];
    T[1] = A[2];
    T[2] = A[1];
    T[3] = A[3];
}

void stranspose_2x2(double *A) {
    const double tmp = A[1];
    A[1] = A[2];
    A[2] = tmp;
}

//! Adds two arrays
void add(const double *x, const double *y, size_t size, double* z) {
    assert( x != NULL && y != NULL && z != NULL );
    for (size_t i = 0; i < size; i++) {
        z[i] = x[i] + y[i];
    }
}

//! Subtracts two arrays
void sub(const double *x, const double *y, size_t size, double* z) {
    assert( x != NULL && y != NULL && z != NULL );
    for (size_t i = 0; i < size; i++) {
        z[i] = x[i] - y[i];
    }
}

//! Scales an array by a scalar
void scal(const double *x, size_t size, double a, double *y) {
    for (size_t i = 0; i < size; i++) {
        y[i] = a*x[i];
    }
}

//! Matrix x Matrix Multiplication: C = A * B 
void mul(const double *A, const double *B, size_t mA, size_t nA, size_t nB, double *C) {
    assert( A != NULL && B != NULL && C != NULL );
    for (size_t i = 0; i < mA; i++) {
        for (size_t j = 0; j < nB; j++) {
            double sum = 0.0;
            for (size_t k = 0; k < nA; k++) {
                sum += A [i*nA+k] * B[k*nB+j];
            } 
            C[i*nB+j] = sum;
        }
    }
}

//! Matrix x Matrix Multiplication ( 2x2 )
void mm_2x2(const double *A, const double *B, double *C) {
    //mm_2x2_avx_v1(A, B, C);
    //mm_2x2_avx_v2(A, B, C);
    C[0] = A[0]*B[0] + A[1]*B[2];
    C[1] = A[0]*B[1] + A[1]*B[3];
    C[2] = A[2]*B[0] + A[3]*B[2];
    C[3] = A[2]*B[1] + A[3]*B[3];
}

//! Matrix x Matrix Multiplication ( 2x2 ) [ AVX ]
void mm_2x2_avx_v1(const double *A, const double *B, double *C) {

    __m256d a = _mm256_load_pd( A );
    __m256d b = _mm256_load_pd( B );
   
    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0101 = _mm256_permute2f128_pd( b, b, 0b00000000 );
    __m256d b2323 = _mm256_permute2f128_pd( b, b, 0b01010101 );

    __m256d c = _mm256_fmadd_pd( a1133, b2323, _mm256_mul_pd( a0022, b0101 ) );
   
    _mm256_store_pd(C, c);
}

//! Matrix x Matrix Multiplication ( 2x2 ) [ AVX ]
void mm_2x2_avx_v2(const double *A, const double *B, double *C) {

    __m256d a = _mm256_load_pd( A );
    __m256d b = _mm256_load_pd( B );

    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0101 = _mm256_permute2f128_pd( b, b, 0b00000000 );
    __m256d b2323 = _mm256_permute2f128_pd( b, b, 0b01010101 );
    
    __m256d c0 = _mm256_mul_pd( a0022, b0101 );
    __m256d c1 = _mm256_mul_pd( a1133, b2323 );

    _mm256_store_pd(C, _mm256_add_pd( c0, c1 ));
}

//! Matrix x Matrix Transpose Multiplication ( 2x2 )
void mmT_2x2(const double *A, const double *B, double *C) {
    //mmT_2x2_avx_v1(A, B, C);
    //mmT_2x2_avx_v2(A, B, C);
    C[0] = A[0]*B[0] + A[1]*B[1];
    C[1] = A[0]*B[2] + A[1]*B[3];
    C[2] = A[2]*B[0] + A[3]*B[1];
    C[3] = A[2]*B[2] + A[3]*B[3];
}

//! Matrix x Matrix Transpose Multiplication ( 2x2 ) [ AVX ]
void mmT_2x2_avx_v1(const double *A, const double *B, double *C) {

    __m256d a = _mm256_load_pd( A );
    __m256d b = _mm256_load_pd( B );

    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0202 = _mm256_permute4x64_pd( b, 0b10001000 );
    __m256d b1313 = _mm256_permute4x64_pd( b, 0b11011101 );
 
    __m256d c = _mm256_fmadd_pd( a1133, b1313, _mm256_mul_pd( a0022, b0202 ) );

    _mm256_store_pd(C, c);
}

//! Matrix x Matrix Transpose Multiplication ( 2x2 ) [ AVX ]
void mmT_2x2_avx_v2(const double *A, const double *B, double *C) {

    __m256d a = _mm256_load_pd( A );
    __m256d b = _mm256_load_pd( B );

    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0202 = _mm256_permute4x64_pd( b, 0b10001000 );
    __m256d b1313 = _mm256_permute4x64_pd( b, 0b11011101 );
    
    __m256d c0 = _mm256_mul_pd( a0022, b0202 );
    __m256d c1 = _mm256_mul_pd( a1133, b1313 );

    _mm256_store_pd(C, _mm256_add_pd( c0, c1 ));
}

//! C += A*B ( 2x2 )
void mmadd_2x2(const double *A, const double *B, double *C) {
    //mmadd_2x2_avx_v1(A, B, C);
    C[0] += A[0]*B[0] + A[1]*B[2];
    C[1] += A[0]*B[1] + A[1]*B[3];
    C[2] += A[2]*B[0] + A[3]*B[2];
    C[3] += A[2]*B[1] + A[3]*B[3];
}

void mmadd_2x2_avx_v1(const double *A, const double *B, double *C) {

    __m256d a = _mm256_load_pd( A );
    __m256d b = _mm256_load_pd( B );
    __m256d c = _mm256_load_pd( C );

    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0101 = _mm256_permute2f128_pd( b, b, 0b00000000 );
    __m256d b2323 = _mm256_permute2f128_pd( b, b, 0b01010101 );
    
    c = _mm256_fmadd_pd( a0022, b0101, c );
    c = _mm256_fmadd_pd( a1133, b2323, c );

    _mm256_store_pd(C, c);
}

void mmadd_2x2_avx_v2(const double *A, const double *B, double *C) {

    __m256d a = _mm256_load_pd( A );
    __m256d b = _mm256_load_pd( B );
    __m256d c = _mm256_load_pd( C );

    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0101 = _mm256_permute2f128_pd( b, b, 0b00000000 );
    __m256d b2323 = _mm256_permute2f128_pd( b, b, 0b01010101 );
    
    c = _mm256_fmadd_pd( a0022, b0101, c );
    
    __m256d cc = _mm256_mul_pd( a1133, b2323 );

    _mm256_store_pd(C, _mm256_add_pd( c, cc ) );
}

//! Matrix x Vector Multiplication ( 2x2 )
void mv_2x2(const double *A, const double *b, double *c) {
    c[0] = A[0]*b[0] + A[1]*b[1];
    c[1] = A[2]*b[0] + A[3]*b[1];
}

//! c += A*b ( 2x2 )
void mvadd_2x2(const double *A, const double *b, double *c) {
    c[0] += A[0]*b[0] + A[1]*b[1];
    c[1] += A[2]*b[0] + A[3]*b[1];
}

//! Cholesky Factorization of a 2x2 SPD Matrix A = L * L^T, L lower triangular
void llt_2x2(const double *A, double *L) {
    assert( A != NULL && L != NULL );
    L[0] = sqrt( A[0] );             // 0*2+0 -> (0,0)
    L[1] = 0.0;                      // 0*2+1 -> (0,1)
    L[2] = A[2] / L[0];              // 1*2+0 -> (1,0)
    L[3] = sqrt( A[3] - L[2]*L[2] ); // 1*2+1 -> (1,1)
}

//! Inverse of a 2x2 Matrix
void inv_2x2(const double *A, double *Ainv) {
    double s = 1.0 / ( A[0]*A[3] - A[1]*A[2] );
    Ainv[0] =  s * A[3];
    Ainv[1] = -s * A[1];
    Ainv[2] = -s * A[2];
    Ainv[3] =  s * A[0];
}

double determinant_2x2(const double* A) {
    return A[0] * A[3] - A[1] * A[2];
}
