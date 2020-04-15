#ifndef TRANSFORMGLOBAL_H
#define TRANSFORMGLOBAL_H
#include "typedefs.h"

/* transform point to global coordinates */
// first rotated and then translated
/*!
    @param[in]  p    2D points
    @param[out] b    x,y translation and rotation in radians
 */
void TransformToGlobal(Matrix2d p, int p_size, Vector3d b);

#endif //TRANSFORMGLOBAL_H
