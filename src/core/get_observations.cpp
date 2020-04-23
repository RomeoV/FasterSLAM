#include "get_observations.h"
#include <iostream>
#include <cmath>

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: TBD
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/

void get_observations(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, int **idf, size_t *nidf, Vector2d (*z))
{
    double *lm_new;
    get_visible_landmarks(x, rmax, lm, lm_rows, &lm_new, idf, nidf); // allocates lm_new
    compute_range_bearing(x, lm_new, *nidf, z);	
    free(lm_new);
}

// lm is a double matrix of dimension nlm x 2
void get_visible_landmarks(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, double **lm_new, int **idf, size_t *nidf)
{
    //select set of landmarks that are visible within vehicle's 
    //semi-circular field of view
    double *dx = (double*) malloc( lm_rows * sizeof(double) );
    double *dy = (double*) malloc( lm_rows * sizeof(double) );

    for (size_t i = 0; i < lm_rows; i++) {
        dx[i] = lm[i*2+0] - x[0];
        dy[i] = lm[i*2+1] - x[1];
    }

    double phi = x[2];

    //distant points are eliminated
    //allocate results of find2() ii and ii_size
    size_t ii_size = lm_rows; 
    size_t *ii = (size_t*) malloc( ii_size * sizeof(size_t) );

    find2(dx, dy, lm_rows, phi, rmax, ii, &ii_size); // fills ii and modifies ii_size

    free(dx); free(dy);

    *lm_new = (double*) malloc( 2*ii_size * sizeof(double) ); // size: 2 x ii_size
    for(size_t i = 0; i < ii_size; i++){
        // lm_new(i,j) = lm(ii[i],j);
        (*lm_new)[ i*2 + 0 ] = lm[ ii[i]*2 + 0 ]; 
        (*lm_new)[ i*2 + 1 ] = lm[ ii[i]*2 + 1 ]; 
    }
 
    int *idf_buff = (int*) malloc( ii_size * sizeof(int) );

    for(int i = 0; i < ii_size; i++) {
        idf_buff[i] = (*idf)[ ii[i] ];
    }

    free(*idf); *idf = idf_buff; idf_buff = NULL; *nidf = ii_size;

    free(ii);
}

// Just fills z which has size lm_rows x 2 and is allocated in the caller of this function
void compute_range_bearing(cVector3d x, const double *lm, const size_t lm_rows, Vector2d (*z)) 
{
    double *dx = (double*) malloc( lm_rows * sizeof(double) ); 
    double *dy = (double*) malloc( lm_rows * sizeof(double) ); 

    for (int i = 0; i < lm_rows; i++) {
        dx[i] = lm[i*2+0] - x[0];
        dy[i] = lm[i*2+1] - x[1];
    }	

    double phi = x[2]; 

    for (int i = 0; i < lm_rows; i++) {
        z[i][0] = sqrt( pow(dx[i],2) + pow(dy[i],2) );
        z[i][1] = atan2( dy[i], dx[i] ) - phi;	
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
