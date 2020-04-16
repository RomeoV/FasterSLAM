#include "data_associate_known.h"

//z is range and bearing of visible landmarks
void data_associate_known(double* z, int* idz, double* table, int Nf, \
						  double* zf, int* idf, double* zn)
{
	// we assume that Nf == number of known landmarks
	// idf.clear();
	for(int i = 0; i < Nf; i++){
		idf[i] = 0;
	}

	int num_observed_features = sizeof(idz) / sizeof(int);
	int idn[Nf]; // dynamic array bound possible?
	int num_new_features = 0;
	int zf_index = 0;
	
	for (int i = 0; i < num_observed_features; i++){ 
		int ii = idz[i]; // check if landmark is in table

		// if this landmark/feature is new --> add to zn
		if (table[ii] == -1) {
			// array is initially empty so we start the index at num_new_features
			zn[num_new_features] = z[i]; // push back to array of known features
			idn[num_new_features] = ii;	
			num_new_features += 1;			
		}
		// unknown --> add to zf
		else {
			// zf array is initially cleared
			zf[zf_index] = z[i]; // push back to array of known features
			int landmark_index = table[ii];
			idf[zf_index] = landmark_index;
		}	
	}

	//assert(idn.size() == zn.size());
	for (int i=0; i < num_new_features; i++) {
		table[idn[i]] = Nf+i;  
	}
}
