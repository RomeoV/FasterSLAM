#pragma once
#include <cstddef>
#include "typedefs.h"

//#ifndef DATA_ASSOCIATE_KNOWN_H
//#define DATA_ASSOCIATE_KNOWN_H

//#include <Eigen/Dense>
//#include <vector>

//using namespace std;
//using namespace Eigen;

/*!
  For each landmark in z, checks in table if this landmark is new (add to zn) or known
  (add to zf). Basically describes data association in continous measurements. [Memory-heavy,
  switch to masks]
  @param[in] 		z 		Landmark observations / measurements [meters, angles].
  @param[in] 		idz		Index of (visible) landmarks.
  @param[out] 	table 	Data association table. -1: New landmark.
  @param[in] 		Nf	 	Number of features / landmarks.
  @param[out] 	zf 		Known landmarks.
  @param[out] 	idf 	Index of known landmarks.
  @param[out] 	zn	 	New landmarks.
  */
void data_associate_known(cVector2d z[], const int* idz, const size_t idz_size, 
        int* table, const int Nf_known, Vector2d zf[], int *idf, size_t *count_zf, Vector2d zn[], size_t *count_zn); 

void data_associate_known_base(cVector2d z[], const int* idz, const size_t idz_size, 
        int* table, const int Nf_known, Vector2d zf[], int *idf, size_t *count_zf, Vector2d zn[], size_t *count_zn);

// Utils
double data_associate_known_base_flops(cVector2d z[], const int* idz, const size_t idz_size, 
                               int* table, const int Nf_known, Vector2d zf[], int *idf, 
                               size_t *count_zf, Vector2d zn[], size_t *count_zn); 

double data_associate_known_active_flops(cVector2d z[], const int* idz, const size_t idz_size, 
                               int* table, const int Nf_known, Vector2d zf[], int *idf, 
                               size_t *count_zf, Vector2d zn[], size_t *count_zn); 

double data_associate_known_base_memory(cVector2d z[], const int* idz, const size_t idz_size, 
                              int* table, const int Nf_known, Vector2d zf[], int *idf, 
                              size_t *count_zf, Vector2d zn[], size_t *count_zn); 

double data_associate_known_active_memory(cVector2d z[], const int* idz, const size_t idz_size, 
                              int* table, const int Nf_known, Vector2d zf[], int *idf, 
                              size_t *count_zf, Vector2d zn[], size_t *count_zn); 
