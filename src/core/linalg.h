#pragma once
#include <iostream>
#include <immintrin.h>

//! ------------------------------------------------------- //
//! ---------- Linear Algebra Utility Functions ----------- //
//! ------------------------------------------------------- //

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

//! Cholesky Factorization of a 2x2 SPD Matrix A = L * L^T, L lower triangular
void llt_2x2(const double *A, double *L);

//! Inverse of a 2x2 Matrix
void inv_2x2(const double *A, double *Ainv);

//! Determinant for 2x2 matrix
//! @param A pointer to row major continuous memory
double determinant_2x2(const double* A);
