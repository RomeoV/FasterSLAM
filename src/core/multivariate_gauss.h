#ifndef MULTIVARIATE_GAUSS_H
#define MULTIVARIATE_GAUSS_H

/*!
    Computes Multivariate Gaussian for n Samples. [Pretty simple adaptions, usually x is size 2 and P is 2x2]
    @param[out]          Sample set. Size len(x) x n (not necessarily a vector if n!=1)
    @param[in]  x        Mean vector (e.g. (V,G)).
    @param[in]  P        Covariance Matrix.
    @param[in]  n        Number of samples.
 */

typedef double* VectorXd;
typedef double* MatrixXd;

void multivariate_gauss(VectorXd x, MatrixXd P, int n, VectorXd result);

#endif //MULTIVARIATE_GAUSS_H
