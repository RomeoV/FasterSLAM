#include "data_associate_known.h"
#include <iostream>

//z is range and bearing of visible landmarks

// vector<Vector2d> z
// vector<int> idz
// VectorXd &table
// int Nf
// vector<Vector2d> &zf
// vector<int> &idf
// vector<Vector2d> &zn
void data_associate_known(const double* z, const int* idz, const size_t idz_size, 
        double* table, const int Nf_known, double *zf, int *idf, size_t *count_zf, double *zn, size_t *count_zn) 
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
            zn[2*(*count_zn)+0] = z[i*2+0];
            zn[2*(*count_zn)+1] = z[i*2+1];
			// idn.push_back(ii);				
            idn[*count_zn] = ii;
            (*count_zn)++;
		}
		else {
			// zf.push_back(z[i]); // z[i] is vector2d
            zf[2*(*count_zf)+0] = z[i*2+0];
            zf[2*(*count_zf)+1] = z[i*2+1];
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
