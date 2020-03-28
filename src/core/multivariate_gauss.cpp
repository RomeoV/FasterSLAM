#include <iostream>
#include "multivariate_gauss.h"
#include "linalg.h"

void multivariate_gauss(VectorXd x, MatrixXd P, int n, VectorXd) 
{
	int len = x.size();
	MatrixXd S = P.llt().matrixL();
	MatrixXd X(len,n);

    fill_rand(X, len*n, -1.0, 1.0);

	MatrixXd ones = MatrixXd::Ones(1,n);	
	return S*X + x*ones;
}
