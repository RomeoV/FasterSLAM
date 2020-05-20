#include "get_observations.h"
#include <math.h>
#include "typedefs.h"
#include "linalg.h"

/*****************************************************************************
 * OPTIMIZATION STATUS
 * Done: Base implementation, unit test
 * ToDo: Start optimizing
 ****************************************************************************/



void get_observations(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, int *idf, size_t *nidf, Vector2d z[]) {
    get_observations_base(x, rmax, lm, lm_rows, idf, nidf, z);
}


/*****************************************************************************
 * PERFORMANCE STATUS
 * Work best: 15*lm_rows flops (get_visible_landmarks)
 * Work worst: 15*lm_rows flops (get_visible_landmarks) + 8*lm_rows flops (range_bearing) = 23 lm_rows
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void get_observations_base(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, int *idf, size_t *nidf, Vector2d z[])
{
    double *lm_new = NULL;
    get_visible_landmarks_base(x, rmax, lm, lm_rows, &lm_new, idf, nidf); // allocates lm_new
    if ( lm_new != NULL ) {
        compute_range_bearing_base(x, lm_new, *nidf, z);	
        free(lm_new);
    }
}


// lm is a double matrix of dimension nlm x 2
void get_visible_landmarks(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, double **lm_new, int *idf, size_t *nidf) {
    get_visible_landmarks_base(x, rmax, lm, lm_rows, lm_new, idf, nidf);
}


/*****************************************************************************
 * PERFORMANCE STATUS
 * Work, best: 13*lm_rows flops (find2) + 2*lm_rows adds = 15*lm_rows flops
 * Work, worst: 13*lm_rows flops (find2) + 2*lm_rows adds = 15*lm_rows flops
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void get_visible_landmarks_base(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, double **lm_new, int *idf, size_t *nidf)
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

    find2_base(dx, dy, lm_rows, phi, rmax, ii, &ii_size); // fills ii and modifies ii_size

    free(dx); free(dy);

    if ( ii_size != 0 ) {

        *lm_new = (double*) malloc( 2*ii_size * sizeof(double) ); // size: 2 x ii_size
        
        for(size_t i = 0; i < ii_size; i++){

            // lm_new(i,j) = lm(ii[i],j);
            (*lm_new)[ i*2 + 0 ] = lm[ ii[i]*2 + 0 ]; 
            (*lm_new)[ i*2 + 1 ] = lm[ ii[i]*2 + 1 ]; 
        }
    
        // use buffer for reordering + copy 
        int *idf_buff = (int*) malloc( ii_size * sizeof(int) );
        // store reordering to buffer
        for(int i = 0; i < ii_size; i++) {
            idf_buff[i] = idf[ ii[i] ];
        }
        // copy back to idf
        for (int i = 0; i < ii_size; i++) {
            idf[i] = idf_buff[i];
        }
        // dealloc buffer 
        free( idf_buff );
    }

    *nidf = ii_size;
    free(ii);
}

void compute_range_bearing(cVector3d x, const double *lm, const size_t lm_rows, Vector2d z[]) {
    compute_range_bearing_base(x, lm, lm_rows, z);
}


/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: 4*lm_rows adds + lm_rows sqrt + lm_rows atan2 + 2*lm_rows pow = 8 flops
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
// Just fills z which has size lm_rows x 2 and is allocated in the caller of this function
void compute_range_bearing_base(cVector3d x, const double *lm, const size_t lm_rows, Vector2d z[]) 
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

void find2(const double *dx, const double *dy, const size_t size, 
        const double phi, const double rmax, size_t *index, size_t *index_size) {
    find2_base(dx,dy,size,phi, rmax, index, index_size);
}

/*****************************************************************************
 * PERFORMANCE STATUS
 * Work, best: 4*size fl-comp + 2*size fabs + 2*size sin + 3*size pow + size adds + 2*size mults = 13 *size flops 
 * Work, worst: 4*size fl-comp + 2*size fabs + 2*size sin + 3*size pow + size adds + 2*size mults = 13 *size flops
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
//! index should be preallocated with size equal to size
void find2_base(const double *dx, const double *dy, const size_t size, 
        const double phi, const double rmax, size_t *index, size_t *index_size)
{
    size_t cnt = 0;
    //incremental tests for bounding semi-circle
    for (size_t j = 0; j < size; j++) {
        const double dxj = dx[j];
        const double dyj = dy[j];
        if ( (fabs(dxj) < rmax) && 
                (fabs(dyj) < rmax) && 
                ((dxj*cos(phi) + dyj*sin(phi)) > 0.0) && 
                ((pow(dxj,2) + pow(dyj,2)) < pow(rmax,2)) )
        {
            index[cnt++] = j;
        }
    }
    *index_size = cnt;
}


//! -----------
//! find2 utils
//! -----------

double find2_base_flops(const double *dx, const double *dy, const size_t size, 
        const double phi, const double rmax, size_t *index, size_t *index_size)
{
    double flops = 0.0;
    // worst case because && is not taken into account
    flops = 4*tp.doublecomp + 2*tp.abs + tp.sin + tp.cos + 2*tp.mul + 2*tp.add + 3*tp.pow;
    return ( flops * size );
}

double find2_active_flops(const double *dx, const double *dy, const size_t size, 
        const double phi, const double rmax, size_t *index, size_t *index_size)
{
    double flops = find2_base_flops(dx, dy, size, phi, rmax, index, index_size);
    return ( flops * size );
}

double find2_base_memory(const double *dx, const double *dy, const size_t size, 
        const double phi, const double rmax, size_t *index, size_t *index_size)
{
    return ( size * 2 );
}

double find2_active_memory(const double *dx, const double *dy, const size_t size, 
        const double phi, const double rmax, size_t *index, size_t *index_size)
{
    double memory = find2_base_memory(dx, dy, size, phi, rmax, index, index_size);
    return memory;
}

//! ---------------------------
//! get_visible_landmarks utils
//! ---------------------------

double get_visible_landmarks_base_flops(cVector3d x, const double rmax, const double *lm, 
        const size_t lm_rows, double **lm_new, int *idf, size_t *nidf)
{
    double flops = 0.0;
    size_t ii_size = 0; // dummy
    flops += lm_rows * 2*tp.add;
    flops += find2_base_flops(NULL, NULL, lm_rows, 0.0, 0.0, NULL, &ii_size);
    return flops;
}

double get_visible_landmarks_active_flops(cVector3d x, const double rmax, const double *lm, 
        const size_t lm_rows, double **lm_new, int *idf, size_t *nidf)
{
    double flops = 0.0;
    size_t ii_size = 0; // dummy
    flops += lm_rows * 2*tp.add;
    flops += find2_active_flops(NULL, NULL, lm_rows, 0.0, 0.0, NULL, &ii_size);
    return flops;
}

double get_visible_landmarks_base_memory(cVector3d x, const double rmax, const double *lm, 
        const size_t lm_rows, double **lm_new, int *idf, size_t *nidf)
{ 
    // ATTENTION!!! *lm_new is not allocated anymore here
    double memory = 0.0;
    
    // select set of landmarks that are visible within vehicle's semi-circular field of view
    double *dx = (double*) malloc( lm_rows * sizeof(double) );
    double *dy = (double*) malloc( lm_rows * sizeof(double) );

    memory += (2*2 + 2*2)*lm_rows + 1;
    for (size_t i = 0; i < lm_rows; i++) {
        dx[i] = lm[i*2+0] - x[0];
        dy[i] = lm[i*2+1] - x[1];
    }

    double phi = x[2];

    // distant points are eliminated, allocate results of find2() ii and ii_size
    size_t ii_size = lm_rows; 
    size_t *ii = (size_t*) malloc( ii_size * sizeof(size_t) );

    find2_base_memory(dx, dy, lm_rows, phi, rmax, ii, &ii_size);
    find2_base(dx, dy, lm_rows, phi, rmax, ii, &ii_size); // fills ii and modifies ii_size

    if ( ii_size != 0 ) {
        memory += ii_size * 15;
    }
    free(dx); free(dy); free(ii); 

    return memory;
}

double get_visible_landmarks_active_memory(cVector3d x, const double rmax, const double *lm, 
        const size_t lm_rows, double **lm_new, int *idf, size_t *nidf)
{ 
    // ATTENTION!!! *lm_new is not allocated anymore here
    return get_visible_landmarks_base_memory(x, rmax, lm, lm_rows, lm_new, idf, nidf);
}

//! ---------------------------
//! compute_range_bearing utils
//! ---------------------------
double compute_range_bearing_base_flops(cVector3d x, const double *lm, const size_t lm_rows, Vector2d z[]) 
{
    return lm_rows * ( 4*tp.add + 2*tp.pow + tp.atan2 + tp.sqrt);
}

double compute_range_bearing_active_flops(cVector3d x, const double *lm, const size_t lm_rows, Vector2d z[]) 
{
    return lm_rows * ( 4*tp.add + 2*tp.pow + tp.atan2 + tp.sqrt);
}

double compute_range_bearing_base_memory(cVector3d x, const double *lm, const size_t lm_rows, Vector2d z[]) 
{
    return 16*lm_rows + 1;
}

double compute_range_bearing_active_memory(cVector3d x, const double *lm, const size_t lm_rows, Vector2d z[]) 
{
    return 16*lm_rows + 1;
}

double get_observations_base_flops(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, int *idf, size_t *nidf, Vector2d z[])
{
    // Save backups for reset
    size_t nidf_bak = *nidf;
    int *idf_bak = (int*) malloc( nidf_bak * sizeof(int) );
    icopy(idf, nidf_bak, idf_bak);

    double flops = 0.0;
    double *lm_new = NULL;

    flops += get_visible_landmarks_base_flops(x, rmax, lm, lm_rows, &lm_new, idf, nidf);

    // explicit call of function to provide correct inputs for compute_range_bearing
    get_visible_landmarks_base(x, rmax, lm, lm_rows, &lm_new, idf, nidf); // allocates lm_new

    if ( lm_new != NULL ) {
        flops += compute_range_bearing_base_flops(x, lm_new, *nidf, z);	
        free(lm_new);
    }
    // Reset 
    icopy(idf_bak, nidf_bak, idf);
    *nidf = nidf_bak;
    free(idf_bak);

    return flops;
}

double get_observations_active_flops(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, int *idf, size_t *nidf, Vector2d z[])
{
    return get_observations_base_flops(x, rmax, lm, lm_rows, idf, nidf, z);
}

double get_observations_base_memory(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, int *idf, size_t *nidf, Vector2d z[])
{
    // Save backups for reset
    size_t nidf_bak = *nidf;
    int *idf_bak = (int*) malloc( nidf_bak * sizeof(int) );
    icopy(idf, nidf_bak, idf_bak);

    double memory = 0.0;
    double *lm_new = NULL;

    memory += get_visible_landmarks_base_memory(x, rmax, lm, lm_rows, &lm_new, idf, nidf);

    // explicit call of function to provide correct inputs for compute_range_bearing
    get_visible_landmarks_base(x, rmax, lm, lm_rows, &lm_new, idf, nidf); // allocates lm_new

    if ( lm_new != NULL ) {
        memory += compute_range_bearing_base_memory(x, lm_new, *nidf, z);	
        free(lm_new);
    }
    // Reset 
    icopy(idf_bak, nidf_bak, idf);
    *nidf = nidf_bak;
    free(idf_bak);

    return memory;
}

double get_observations_active_memory(cVector3d x, const double rmax, const double *lm, const size_t lm_rows, int *idf, size_t *nidf, Vector2d z[])
{
    return get_observations_base_memory(x, rmax, lm, lm_rows, idf, nidf, z);
}
