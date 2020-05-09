#include "linalg.h"
#include <iostream>
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
    C[0] = A[0]*B[0] + A[1]*B[2];
    C[1] = A[0]*B[1] + A[1]*B[3];
    C[2] = A[2]*B[0] + A[3]*B[2];
    C[3] = A[2]*B[1] + A[3]*B[3];
}

//! Matrix x Vector Multiplication ( 2x2 )
void mv_2x2(const double *A, const double *B, double *C) {
    C[0] = A[0]*B[0] + A[1]*B[1];
    C[1] = A[2]*B[0] + A[3]*B[1];
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
