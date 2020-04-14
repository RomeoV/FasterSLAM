#include "get_observations.h"
#include <iostream>
#include <cmath>

void get_observations(cVector3d x, double *lm, size_t lm_cols, int **idf, size_t *nidf, const double rmax, double *z)
{
	get_visible_landmarks(x, &lm, &lm_cols, idf, nidf, rmax);
	compute_range_bearing(x, lm, lm_cols, z);	
}

// lm is a double matrix of dimension 2 x nlm
// !!! Pointers and sizes of ``lm'' and ``idf'' are modified inside get_visible_landmarks() !!!
void get_visible_landmarks(cVector3d x, double **lm, size_t *lm_cols, int **idf, size_t *nidf, const double rmax)
{
	//select set of landmarks that are visible within vehicle's 
	//semi-circular field of view
	double *dx = (double*) malloc( *lm_cols * sizeof(double) );
	double *dy = (double*) malloc( *lm_cols * sizeof(double) );

	for (size_t j = 0; j < *lm_cols; j++) {
		dx[j] = (*lm)[0*(*lm_cols)+j] - x[0];
		dy[j] = (*lm)[1*(*lm_cols)+j] - x[1];
	}

	double phi = x[2];

	//distant points are eliminated
    //allocate results of find2() ii and ii_size
    size_t ii_size = *lm_cols; 
    size_t *ii = (size_t*) malloc( ii_size * sizeof(size_t) );

    find2(dx, dy, *lm_cols, phi, rmax, ii, &ii_size); // fills ii and modifies ii_size

    free(dx); free(dy);

	double *lm_new = (double*) malloc( 2*ii_size * sizeof(double) ); // size: 2 x ii_size
	for (size_t j = 0; j < 2; j++){
		for(size_t k = 0; k < ii_size; k++){
            // lm_new(j,k) = lm(j,ii[k]);
            lm_new[ j*ii_size + k ] = (*lm)[ j*(*lm_cols) + ii[k] ]; 
		}
	}

	// Do lm = MatrixXd(lm_new) without reallocation
    free(*lm); *lm = lm_new; lm_new = NULL; *lm_cols = ii_size;

    int *idf_buff = (int*) malloc( ii_size * sizeof(int) );

	for(int i=0; i<ii_size; i++) {
		idf_buff[i] = (*idf)[ ii[i] ];
	}

    free(*idf); *idf = idf_buff; idf_buff = NULL; *nidf = ii_size;
    free(ii);
}

// Just fills z which has size lm_cols x 2 and is allocated in the caller of this function
void compute_range_bearing(cVector3d x, const double *lm, const size_t lm_cols, double *z) 
{
	double *dx = (double*) malloc( lm_cols * sizeof(double) ); 
	double *dy = (double*) malloc( lm_cols * sizeof(double) ); 

	for (int j = 0; j < lm_cols; j++) {
		dx[j] = lm[0*lm_cols+j] - x[0];
		dy[j] = lm[1*lm_cols+j] - x[1];
	}	

	double phi = x[2]; 

	for (int j = 0; j < lm_cols; j++) {
		z[j*2 + 0] = sqrt( pow(dx[j],2) + pow(dy[j],2) );
        z[j*2 + 1] = atan2( dy[j], dx[j] ) - phi;	
	}

    free(dx); free(dy);
}

//! index should be preallocated with size equal to size
void find2(const double *dx, const double *dy, const size_t size, 
           const double phi, const double rmax, size_t *index, size_t *index_size)
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
