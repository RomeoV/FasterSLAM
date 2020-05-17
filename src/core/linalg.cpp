#include "linalg.h"
#include <iostream>
#include <immintrin.h>
#include <cmath>
#include <cassert>
#include <immintrin.h>
#include "fastrand.h"

#include <stdint.h>
#include <string.h>

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

void print256d(__m256d var)
{
    double val[4];
    memcpy(val, &var, sizeof(val));
    printf("SIMD: <%f %f %f %f> \n", 
           val[0], val[1], val[2], val[3]);
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
    // assert( x != NULL && hi > lo );
    double range = hi - lo;
    for (size_t i = 0; i < size; i++) {
        x[i] = lo + range * ((double)rand())/((double)RAND_MAX);
    }
}



// Works ;)
#ifdef __AVX2__
__m256d fill_rand_avx(double lo, double hi) {
    __m256d lov = _mm256_set1_pd(lo);
    __m256d hiv = _mm256_set1_pd(hi);
    __m256i rand_vec =  avx_xorshift128plus();
    
    __m256d intmax = _mm256_set1_pd(2147483647);
    __m256d uintmax = _mm256_set1_pd(4294967295);

    // We only use 128bit of the vector atm. E.g. for predict-update, we could make use of the other 128 bit too.
    __m256d ymm0 = _mm256_cvtepi32_pd(_mm256_extractf128_si256(rand_vec, 0)); //32 bit
    ymm0 = _mm256_add_pd(ymm0, intmax);
    __m256d ymm1 = _mm256_div_pd(ymm0, uintmax);

    __m256d range = _mm256_sub_pd(hiv, lov);
#ifdef __FMA__
    return _mm256_fmadd_pd(range, ymm1, lov);
#else
    return _mm256_add_pd( _mm256_mul_pd(range, ymm1), lov);
#endif
}
#endif

// Fails some tests
#ifdef __AVX2__
__m256d fill_rand_avx_abs(double lo, double hi) {
    __m256d lov = _mm256_set1_pd(lo);
    __m256d hiv = _mm256_set1_pd(hi);
    __m256i rand_vec =  _mm256_abs_epi32(avx_xorshift128plus());
    
    __m256d intmax = _mm256_set1_pd(2147483647); //2147483647
    __m256d ymm0 = _mm256_cvtepi32_pd(_mm256_extractf128_si256(rand_vec, 0)); //32 bit
    __m256d ymm1 = _mm256_div_pd(ymm0, intmax);

    __m256d range = _mm256_sub_pd(hiv, lov);
#ifdef __FMA__
    return _mm256_fmadd_pd(range, ymm1, lov);
#else
    return _mm256_add_pd( _mm256_mul_pd(range, ymm1), lov);
#endif
}
#endif

/*
void fill_rand_fast(double *x, size_t size, double lo, double hi) {
    double range = hi - lo;
    for (size_t i = 0; i < size; i++) {
        x[i] = lo + range * xorshf96()*1.0/ulong_max;
    }
}
*/

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
#ifdef __AVX__
void mm_2x2_avx_v1(const double *A, const double *B, double *C) {

    __m256d a = _mm256_load_pd( A );
    __m256d b = _mm256_load_pd( B );
   
    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0101 = _mm256_permute2f128_pd( b, b, 0b00000000 );
    __m256d b2323 = _mm256_permute2f128_pd( b, b, 0b01010101 );

#ifdef __FMA__
    __m256d c = _mm256_fmadd_pd( a1133, b2323, _mm256_mul_pd( a0022, b0101 ) );
#else
    __m256d c = _mm256_add_pd( _mm256_mul_pd(a1133, b2323), _mm256_mul_pd( a0022, b0101 ) );
#endif
   
    _mm256_store_pd(C, c);
}
#endif

//! Matrix x Matrix Multiplication ( 2x2 ) [ AVX ]
#ifdef __AVX__
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
#endif

#ifdef __AVX__
void mm_2x2_avx_v3(const double *A, const double *B, double *C) {
     
    __m256d a0123 = _mm256_load_pd( A );
    __m256d b0123 = _mm256_load_pd( B );
    
    __m256d a1032 = _mm256_permute_pd( a0123, 0b0101 );
    __m256d b2301 = _mm256_permute2f128_pd( b0123, b0123, 0b00000001 );
    __m256d b0303 = _mm256_blend_pd( b0123, b2301, 0b0110 );
    __m256d b2121 = _mm256_blend_pd( b0123, b2301, 0b1001 );

#ifdef __FMA__
    __m256d c = _mm256_fmadd_pd( a1032, b2121, _mm256_mul_pd( a0123, b0303 ) );
#else
    __m256d c = _mm256_add_pd( _mm256_mul_pd(a1032, b2121), _mm256_mul_pd( a0123, b0303 ) );
#endif

    _mm256_store_pd(C, c);
}
#endif

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
#ifdef __AVX2__
void mmT_2x2_avx_v1(const double *A, const double *B, double *C) {

    __m256d a = _mm256_load_pd( A );
    __m256d b = _mm256_load_pd( B );

    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0202 = _mm256_permute4x64_pd( b, 0b10001000 );
    __m256d b1313 = _mm256_permute4x64_pd( b, 0b11011101 );
 
#ifdef __FMA__
    __m256d c = _mm256_fmadd_pd( a1133, b1313, _mm256_mul_pd( a0022, b0202 ) );
#else
    __m256d c = _mm256_add_pd( _mm256_mul_pd(a1133, b1313), _mm256_mul_pd( a0022, b0202 ) );
#endif

    _mm256_store_pd(C, c);
}
#endif

//! Matrix x Matrix Transpose Multiplication ( 2x2 ) [ AVX ]
#ifdef __AVX2__
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
#endif

#ifdef __AVX__
void mmT_2x2_avx_v3(const double *A, const double *B, double *C) {

    __m256d a0123 = _mm256_load_pd( A );
    __m256d b0123 = _mm256_load_pd( B );
    
    __m256d b2301 = _mm256_permute2f128_pd( b0123, b0123, 0b00000001 );
    __m256d b0303 = _mm256_blend_pd( b0123, b2301, 0b0110 );
    __m256d b2121 = _mm256_blend_pd( b0123, b2301, 0b1001 );
    
    __m256d c0 = _mm256_mul_pd( a0123, b2121 );
    __m256d c0_perm = _mm256_permute_pd( c0, 0b0101 );
#ifdef __FMA__
    __m256d c = _mm256_fmadd_pd( a0123, b0303, c0_perm );
#else
    __m256d c = _mm256_add_pd( _mm256_mul_pd(a0123, b0303), c0_perm );
#endif

    _mm256_store_pd( C, c );
}
#endif



//! C += A*B ( 2x2 )
void mmadd_2x2(const double *A, const double *B, double *C) {
    //mmadd_2x2_avx_v1(A, B, C);
    C[0] += A[0]*B[0] + A[1]*B[2];
    C[1] += A[0]*B[1] + A[1]*B[3];
    C[2] += A[2]*B[0] + A[3]*B[2];
    C[3] += A[2]*B[1] + A[3]*B[3];
}

#ifdef __AVX__
void mmadd_2x2_avx_v1(const double *A, const double *B, double *C) {

    __m256d a = _mm256_load_pd( A );
    __m256d b = _mm256_load_pd( B );
    __m256d c = _mm256_load_pd( C );

    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0101 = _mm256_permute2f128_pd( b, b, 0b00000000 );
    __m256d b2323 = _mm256_permute2f128_pd( b, b, 0b01010101 );
    
#ifdef __FMA__
    c = _mm256_fmadd_pd( a0022, b0101, c );
#else
    c = _mm256_add_pd( _mm256_mul_pd(a0022, b0101), c );
#endif
#ifdef __FMA__
    c = _mm256_fmadd_pd( a1133, b2323, c );
#else
    c = _mm256_add_pd( _mm256_mul_pd(a1133, b2323), c );
#endif

    _mm256_store_pd(C, c);
}
#endif

#ifdef __AVX__
void mmadd_2x2_avx_v2(const double *A, const double *B, double *C) {

    __m256d a = _mm256_load_pd( A );
    __m256d b = _mm256_load_pd( B );
    __m256d c = _mm256_load_pd( C );

    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0101 = _mm256_permute2f128_pd( b, b, 0b00000000 );
    __m256d b2323 = _mm256_permute2f128_pd( b, b, 0b01010101 );
    
#ifdef __FMA__
    c = _mm256_fmadd_pd( a0022, b0101, c );
#else
    c = _mm256_add_pd( _mm256_mul_pd(a0022, b0101), c );
#endif
    
    __m256d cc = _mm256_mul_pd( a1133, b2323 );

    _mm256_store_pd(C, _mm256_add_pd( c, cc ) );
}
#endif

#ifdef __AVX__
__m256d _mmadd_2x2_avx_v2(__m256d a, __m256d b, __m256d c) {

    __m256d a0022 = _mm256_permute_pd( a, 0b0000 );
    __m256d a1133 = _mm256_permute_pd( a, 0b1111 );
    __m256d b0101 = _mm256_permute2f128_pd( b, b, 0b00000000 );
    __m256d b2323 = _mm256_permute2f128_pd( b, b, 0b01010101 );
    
#ifdef __FMA__
    c = _mm256_fmadd_pd( a0022, b0101, c );
#else
    c = _mm256_add_pd( _mm256_mul_pd(a0022, b0101), c );
#endif
    
    __m256d cc = _mm256_mul_pd( a1133, b2323 );

    return _mm256_add_pd( c, cc );
}
#endif

#ifdef __AVX__
__m256d _mmTadd_2x2_avx_v2(__m256d a0123, __m256d b0123, __m256d c) {

    __m256d b2301 = _mm256_permute2f128_pd( b0123, b0123, 0b00000001 );
    __m256d b0303 = _mm256_blend_pd( b0123, b2301, 0b0110 );
    __m256d b2121 = _mm256_blend_pd( b0123, b2301, 0b1001 );
    
    __m256d c0 = _mm256_mul_pd( a0123, b2121 );
#ifdef __FMA__
    __m256d cc = _mm256_fmadd_pd( a0123, b0303, c );
#else
    __m256d cc = _mm256_add_pd( _mm256_mul_pd(a0123, b0303), c );
#endif
    __m256d c0_perm = _mm256_permute_pd( c0, 0b0101 );

    return _mm256_add_pd(cc, c0_perm);
}
#endif

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

#ifdef __AVX2__
/** v.T @ M @ v subroutine
 * @param m1 2x2 Matrix M1 in row major storage
 * @param m2 2x2 Matrix M2 in row major storage
 * @param v12 The 2x1 vectors V1 and V2 stored in order
 *
 * @returns AVX2[a, b , c, d] s.t. a + b = V1.T@M1@V1, c + d = V2.T@M2@V2
 */
__m256d produce_hvec(const __m256d& m1, const __m256d& m2, const __m256d& v12) {
  __m256d v12_v12 = _mm256_mul_pd(v12, v12);  // [v1^2, v2^2, v3^2, v4^2]
  __m256d v12_v21 = _mm256_mul_pd(v12, _mm256_permute_pd(v12, 0b0101));  // [v1v2, v1v2, v3v4, v3v4]

  __m256d m1_perm = _mm256_permute4x64_pd(m1, 0b10011100);  // [a, d, b, c]
  __m256d m2_perm = _mm256_permute4x64_pd(m2, 0b11001001);  // [e, h, f, g
  __m256d m1_     = _mm256_permute2f128_pd(m1_perm, m2_perm, 0b00110000);  // [a, d, e, h]
  __m256d m2_     = _mm256_permute2f128_pd(m1_perm, m2_perm, 0b00100001);  // [b, c, f, g]

  __m256d lhs      = _mm256_mul_pd(v12_v12, m1_);  // [a11, d22, e33, h44]
  __m256d rhs      = _mm256_mul_pd(v12_v21, m2_);  // [b12, c12, f34, g34]
  __m256d res      = _mm256_add_pd(lhs, rhs);  // [a11+b12, c12+d22, e33+f34, g34+h44]
                                               // ^ needs to be hsum'd to get final result

  return res;
};
#endif

#ifdef __AVX2__
__m256d mm_vT_M_v_avx2(const __m256d& m1,  const __m256d& m2,
                       const __m256d& m3,  const __m256d& m4,
                       const __m256d& v12, const __m256d& v34) {
  __m256d res1 = produce_hvec(m1, m2, v12);
  __m256d res2 = produce_hvec(m3, m4, v34);

  __m256d res_total = _mm256_hadd_pd(res1, res2);  // this gives [A, C, B, D] so we have to permute B and C
  __m256d result = _mm256_permute4x64_pd(res_total, 0b11011000);
  return result;
}
#endif
