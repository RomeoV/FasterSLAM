#pragma once
#include <cstddef>

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
//void data_associate_known(vector<Vector2d> z, vector<int> idz, VectorXd &table, int Nf, \
//						  vector<Vector2d> &zf, vector<int> &idf, vector<Vector2d> &zn); 

void data_associate_known(const double* z, const int* idz, const size_t idz_size, 
        double* table, const int Nf_known, double *zf, int *idf, size_t *count_zf, double *zn, size_t *count_zn); 
