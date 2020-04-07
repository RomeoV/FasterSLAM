#include "get_observations.h"
#include <iostream>
#include <cmath>

vector<Vector2d> get_observations(Vector3d x, MatrixXd lm, vector<int> &idf, double rmax)
{
	get_visible_landmarks(x, lm, idf, rmax);
	return compute_range_bearing(x, lm);	
}


void get_visible_landmarks(Vector3d x, MatrixXd &lm, vector<int> &idf, double rmax)
{
	//select set of landmarks that are visible within vehicle's 
	//semi-circular field of view
	vector<double> dx;
	vector<double> dy;

	for (int c=0; c<lm.cols(); c++) {
		dx.push_back(lm(0,c) - x(0));
		dy.push_back(lm(1,c) - x(1));
	}

	double phi = x(2);

	//distant points are eliminated
	vector<int> ii = find2(dx,dy,phi,rmax); 

	MatrixXd lm_new (lm.rows(), ii.size());
	unsigned j,k;
	for (j=0; j<lm.rows(); j++){
		for(k=0; k< ii.size(); k++){
			lm_new(j,k) = lm(j,ii[k]);
		}
	}
	lm = MatrixXd(lm_new); 

	vector<int> idf_backup(idf);
	idf.clear();
	for(int i=0; i<ii.size(); i++) {
		idf.push_back(idf_backup[ii[i]]);
	}
}

vector<Vector2d> compute_range_bearing(Vector3d x, MatrixXd lm) 
{
	vector<double> dx; 
	vector<double> dy;

	for (int c=0; c<lm.cols(); c++) {
		dx.push_back(lm(0,c) - x(0));
		dy.push_back(lm(1,c) - x(1));
	}	

	assert(dx.size() == lm.cols());
	assert(dy.size() == lm.cols());

	double phi = x(2);
	vector<Vector2d> z;

	for (int i =0; i<lm.cols(); i++) {
		Vector2d zvec(2);
		zvec<< sqrt(pow(dx[i],2) + pow(dy[i],2)), atan2(dy[i],dx[i]) - phi;	
		z.push_back(zvec);
	}

	return z; 
}
//! index should be preallocated with size equal to size
void find2(const double *dx, const double *dy, const size_t size, 
           const double phi, const double rmax, double *index, size_t *index_size)
{
    size_t cnt = 0;
	//incremental tests for bounding semi-circle
	for (size_t j = 0; j < size; j++) {
        const double dxj = dx[j];
        const double dyj = dy[j];
		if ( (abs(dxj) < rmax) && 
             (abs(dyj) < rmax) && 
             ((dxj*cos(phi) + dyj*sin(phi)) > 0.0) && 
             ((pow(dxj,2) + pow(dyj,2)) < pow(rmax,2)) )
        {
            index[cnt++] = j;
        }
	}
    *index_size = cnt;
}
