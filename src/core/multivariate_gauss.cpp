#include <iostream>
#include "multivariate_gauss.h"
#include "linalg.h"

void multivariate_gauss(VectorXd x, MatrixXd P, int n, VectorXd) 
{
	int len = x.size();
	MatrixXd S = P.llt().matrixL();
	MatrixXd X(len,n);
	
    double LO = -1.0f;
    double HI = 1.0f;

    for (int i = 0; i < len; i++) {
        for (int j = 0; j < n; j++) {
            double r3 = LO + (double)rand()/((double)RAND_MAX/(HI-LO));
            X(i,j) = r3;
        }
    }
	
	MatrixXd ones = MatrixXd::Ones(1,n);	
	return S*X + x*ones;
}
