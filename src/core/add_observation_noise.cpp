#include "add_observation_noise.h"
#if 0
//http://moby.ihme.washington.edu/bradbell/mat2cpp/randn.cpp.xml
MatrixXd randn(int m, int n) 
{	
	// use formula 30.3 of Statistical Distributions (3rd ed)
	// Merran Evans, Nicholas Hastings, and Brian Peacock
	int urows = m * n + 1;
	MatrixXd u(urows, 1);

	//u is a random matrix
	for (int r=0; r<urows; r++) {
		for (int c=0; c<1; c++) {
			u(r,c) = std::rand();
		}
	}
	
	MatrixXd x(m,n);

	int i, j, k;
	double pi = 4. * std::atan(1.f);
	double square, amp, angle;
	k = 0;
	for(i = 0; i < m; i++)
	{	for(j = 0; j < n; j++)
		{	if( k % 2 == 0 )
			{	square = - 2. * std::log( u(k, 0) );
				if( square < 0. )
					square = 0.;
				amp = sqrt(square);
				angle = 2. * pi * u(k+1, 0);
				x(i, j) = amp * std::sin( angle );
			}
			else	
				x(i, j) = amp * std::cos( angle );
			k++;
		}
	}
	return x;
}

MatrixXd rand(int m, int n) 
{	
	MatrixXd x(m,n);	
	int i, j;
	double rand_max = double(RAND_MAX);

	for(i = 0; i < m; i++) {	
		for(j = 0; j < n; j++)
			x(i, j) = double(std::rand()) / rand_max;
	} 
	return x;
}
#endif
//add random measurement noise. We assume R is diagnoal matrix
void add_observation_noise(vector<Vector2d> &z, Matrix2d R, int addnoise)
{
    double LO = -1.0f;
    double HI = 1.0f;

	if (addnoise == 1){
		int len = z.size();	
		if (len > 0) {
			//MatrixXd randM1 = nRandMat::randn(1,len);
            MatrixXd randM1(1,len);
            for (int i=0; i< len; i++) {
                double r3 = LO + (double)rand()/((double)RAND_MAX/(HI-LO));
                randM1(0,i) = r3; 
            }
			//MatrixXd randM2 = nRandMat::randn(1,len);
			MatrixXd randM2(1,len);
            for (int j=0; j< len; j++) {
                double r4 = LO + (double)rand()/((double)RAND_MAX/(HI-LO));
                randM2(0,j) = r4; 
            }
            //cout<<"randM1"<<endl;
            //cout<<randM1<<endl;

			for (int c=0; c<len; c++) {
				z[c][0] = z[c][0] + randM1(0,c)*sqrt(R(0,0));
				z[c][1] = z[c][1] + randM2(0,c)*sqrt(R(1,1));
			}
		}
	}	
}


