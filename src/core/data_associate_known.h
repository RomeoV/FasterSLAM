#pragma once
#include "typedefs.h"

/*!
	For each landmark in z, checks in table if this landmark is new (add to zn) or known
	(add to zf). Basically describes data association in continous measurements. [Memory-heavy,
	switch to masks]
	@param[in] 		z 		Landmark observations / measurements [meters, angles].
	@param[in] 		idz		Index of (visible) landmarks.
	@param[out] 	table 	Data association table. -1: New landmark.
	@param[in] 		Nf	 	Number of features / landmarks. /
		--> seems this is the number of features in the table/known landmarks
	@param[out] 	zf 		Known landmarks.
	@param[out] 	idf 	Index of known landmarks.
	@param[out] 	zn	 	New landmarks.
 */
// z is vector<Vector2d>
// idz is vector<int>
// &table is VectorXd
void data_associate_known(double* z, int* idz, double* table, int Nf, 
                          double* zf, int* idf, double* zn);
