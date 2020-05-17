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


// Source: https://latkin.org/blog/2014/11/09/a-simple-benchmark-of-various-math-operations/, Intel Intrinsics guide
typedef struct throughputs_s {
    double sin = 15.0; //All values set to 1.0 for now, we change later
    double atan2 = 25.0;
    double mul = 1.0;
    double add = 1.0;
    double div = 5.0;
    double sqrt = 7.0;
    double rsqrt = 1.0;
    double negation = 1.0;
    double cos = 15.0;
    double doublecomp = 1.0;
    double rand = 30.0; //Guess
    double fastrand = 7.0; //Guess
    double floor = 1.0;
    double modulo = 15.0; //Guess
    double pow = 1.0; // 1 Mult  = pow2 (The only one we use)
    double exp = 10.0;
} throughputs_t;

static throughputs_t tp;

