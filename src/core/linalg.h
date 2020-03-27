#ifndef LINALG_H
#define LINALG_H

void fill(double *x, int size, double val) {
    for (int i = 0; i < size; i++) {
        x[i] = val; 
    }
}

void fill_rand(double *x, int size, int lo, int hi) {
    for (int i = 0; i < size; i++) {
    
    }
}

void add(double *x, double *y, int size, double* z) {
    for (int i = 0; i < size; i++) {
        z[i] = x[i] + y[i];
    }
}

// Matrix x Matrix Multiplication: C = A * B 
// Matrix A ( mA x nA )
// Matrix B ( nA x nB )
// Matrix C ( mA x nB )
void mul(double *A, double *B, int mA, int nA, int nB, double *C) {
    for (int i = 0; i < mA; i++) {
        for (int j = 0; j < nB; j++) {
            double sum = 0.0;
            for (int k = 0; k < nA; k++) {
                sum += A [i*nA+k] * B[k*nB+j];
            } 
            C[i*nB+j] = sum;
        }
    }
}

// Cholesky Factorization A = L * L^T
// Matrix A ( mA x mA )
// Matrix L ( mA x mA, Lower triangular ) 
void llt(double *A, int mA, double *L) {

}

#endif LINALG_H
