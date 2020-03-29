#pragma once
#include <iostream>

//! ------------------------------------------------------- //
//! ---------- Linear Algebra Utility Functions ----------- //
//! ------------------------------------------------------- //

//! Prints an array
void print(double *x, size_t rows, size_t cols);

//! Fills an array with a specific value
void fill(double *x, size_t size, double val);

//! Fills an array with random values in the range [lo, hi]
void fill_rand(double *x, size_t size, double lo, double hi);

//! ------------------------------------------------------- //
//! -------------- Basic Matrix Operations ---------------- //
//! ------------------------------------------------------- //

//! Matrix Transpose
void transpose(double *A, size_t mA, size_t nA, double *T); 

//! Adds two arrays
void add(double *x, double *y, size_t size, double* z); 

//! Subtracts two arrays
void sub(double *x, double *y, size_t size, double* z);

//! Scales an array by a scalar
void scal(double *x, size_t size, double a);

//! Matrix x Matrix Multiplication: C = A * B 
void mul(const double *A, const double *B, size_t mA, size_t nA, size_t nB, double *C);

//! Cholesky Factorization of a 2x2 SPD Matrix A = L * L^T, L lower triangular
void llt_2x2(double *A, double *L);

//! Inverse of a 2x2 Matrix
void inv_2x2(const double *A, double *Ainv);
