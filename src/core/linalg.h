// #include "linalg.h"
#pragma once
#include <iostream>
#include <cassert>

// ------------------------------------------------------- //
// ---------- Linear Algebra Utility Functions ----------- //
// ------------------------------------------------------- //

// Prints an array
void print(double *x, size_t rows, size_t cols) {
    assert( x != NULL );
    for (size_t i = 0; i < rows; i++) {
        std::cout << std::endl;
        for (size_t j = 0; j < cols; j++) {
            std::cout << "\t" << x[i*cols+j];
        }
    }
    std::cout << std::endl;
}

// Fills an array with a specific value
void fill(double *x, size_t size, double val) {
    assert( x != NULL );
    for (size_t i = 0; i < size; i++) {
        x[i] = val; 
    }
}

// Fills an array with random values in the range [lo, hi]
void fill_rand(double *x, size_t size, double lo, double hi) {
    assert( x != NULL && hi > lo );
    double range = hi - lo;
    for (size_t i = 0; i < size; i++) {
        x[i] = lo + range * ((double)rand())/((double)RAND_MAX);
    }
}

// ------------------------------------------------------- //
// -------------- Basic Matrix Operations ---------------- //
// ------------------------------------------------------- //

// Matrix Transpose
// Matrix A ( mA x nA ) -> Input
// Matrix T ( nA x mA ) -> Output
// CAUTION!!! For in place transposition use stranspose() instead
void transpose(double *A, size_t mA, size_t nA, double *T) {
    assert( A != NULL && T != NULL );
    for (size_t i = 0; i < nA; i++) {
        for (size_t j = 0; j < mA; j++) {
             T[i*mA+j] = A[j*nA+i];
        }
    }
}

// Adds two arrays
void add(double *x, double *y, size_t size, double* z) {
    assert( x != NULL && y != NULL && z != NULL );
    for (size_t i = 0; i < size; i++) {
        z[i] = x[i] + y[i];
    }
}

// Subtracts two arrays
void sub(double *x, double *y, size_t size, double* z) {
    assert( x != NULL && y != NULL && z != NULL );
    for (size_t i = 0; i < size; i++) {
        z[i] = x[i] - y[i];
    }
}

// Matrix x Matrix Multiplication: C = A * B 
// Matrix A ( mA x nA ) -> Input
// Matrix B ( nA x nB ) -> Input
// Matrix C ( mA x nB ) -> Output
void mul(double *A, double *B, size_t mA, size_t nA, size_t nB, double *C) {
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

// Cholesky Factorization A = L * L^T
// Matrix A ( mA x mA )
// Matrix L ( mA x mA, Lower triangular ) 
void llt(double *A, size_t mA, double *L) {
    assert( A != NULL && L != NULL );
}
