#pragma once
#include <cstdlib>

typedef double Vector2d[2];
typedef double Vector3d[3];
typedef double Matrix2d[4];
typedef double Matrix3d[9];
typedef double Matrix23d[6];
typedef const double cVector2d[2];
typedef const double cVector3d[3];
typedef const double cMatrix2d[4];
typedef const double cMatrix3d[9];
typedef const double cMatrix23d[6];

typedef struct throughputs_s {
    double sin = 1.0; //All values set to 1.0 for now, we change later
    double atan2 = 1.0;
    double mul = 1.0;
    double add = 1.0;
    double div = 1.0;
    double sqrt = 1.0;
    double rsqrt = 1.0;
    double negation = 1.0;
    double cos = 1.0;
    double doublecomp = 1.0;
} throughputs_t;

