#include "data_associate_known.h"
#include <iostream>
#include "typedefs.h"

//z is range and bearing of visible landmarks

// vector<Vector2d> z
// vector<int> idz
// VectorXd &table
// int Nf
// vector<Vector2d> &zf
// vector<int> &idf
// vector<Vector2d> &zn

void data_associate_known(cVector2d z[], const int* idz, const size_t idz_size, 
        int* table, const int Nf_known, Vector2d zf[], int *idf, size_t *count_zf, Vector2d zn[], size_t *count_zn) {
    data_associate_known_base(z, idz, idz_size, table,Nf_known, zf, idf, count_zf, zn, count_zn);
}



/*****************************************************************************
 * PERFORMANCE STATUS
 * Work: 0 flops + count_zn int-adds = 0 flops
 * Memory moved: TBD
 * Cycles: TBD
 * Performance: TBD
 * Optimal: TBD
 * Status: TBD
 ****************************************************************************/
void data_associate_known_base(cVector2d z[], const int* idz, const size_t idz_size, 
        int* table, const int Nf_known, Vector2d zf[], int *idf, size_t *count_zf, Vector2d zn[], size_t *count_zn) 
{
    // zn and zf are always allocated but considered empty in this step
    // idf.clear(); // dealloc or just set to zero 
    int *idn = (int*) malloc(idz_size * sizeof(int));

    int i = 0, ii = 0;
    *count_zn = 0;
    *count_zf = 0; 
    for (i = 0; i < idz_size; i++){
        ii = idz[i];
        if ( table[ii] == -1 ) { // new feature
            // zn.push_back(z[i]); // z[i] is vector2d
            zn[*count_zn][0] = z[i][0];
            zn[*count_zn][1] = z[i][1];
            // idn.push_back(ii);				
            idn[*count_zn] = ii;
            (*count_zn)++;
        }
        else {
            // zf.push_back(z[i]); // z[i] is vector2d
            zf[*count_zf][0] = z[i][0];
            zf[*count_zf][1] = z[i][1];
            // idf.push_back( table[ii] );
            idf[*count_zf] = table[ii];
            (*count_zf)++;
        }	
    }

    // assert(idn.size() == zn.size());
    for (int i = 0; i < *count_zn; i++) {
        table[ idn[i] ] = Nf_known + i;  
    }
    free(idn);
}

// Utils

FlopCount data_associate_known_base_flops(cVector2d z[], const int* idz, const size_t idz_size, 
                               int* table, const int Nf_known, Vector2d zf[], int *idf, 
                               size_t *count_zf, Vector2d zn[], size_t *count_zn) 
{
    return FlopCount();  // it's literally zero
}

FlopCount data_associate_known_active_flops(cVector2d z[], const int* idz, const size_t idz_size, 
                               int* table, const int Nf_known, Vector2d zf[], int *idf, 
                               size_t *count_zf, Vector2d zn[], size_t *count_zn) 
{
    return FlopCount();  // it's literally zero
}

double data_associate_known_base_memory(cVector2d z[], const int* idz, const size_t idz_size, 
                              int* table, const int Nf_known, Vector2d zf[], int *idf, 
                              size_t *count_zf, Vector2d zn[], size_t *count_zn) 
{
    return ( 2*2 + 2 ) * idz_size;
}

double data_associate_known_active_memory(cVector2d z[], const int* idz, const size_t idz_size, 
                              int* table, const int Nf_known, Vector2d zf[], int *idf, 
                              size_t *count_zf, Vector2d zn[], size_t *count_zn) 
{
    return ( 2*2 + 2 ) * idz_size;
}
