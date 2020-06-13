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
// See flop_estimates.cpp for a benchmark of these realtive costs, try different input scales there
typedef struct throughputs_s {
    double abs = 4.0;
    double sin = 19.0;
    double atan2 = 30.0;
    double mul = 1.0;
    double add = 1.0;
    double div = 5.0;
    double sqrt = 7.0;
    double rsqrt = 1.0; // Not used atm
    double negation = 1.0;
    double cos = 21.0;
    double doublecomp = 10.0;
    double rand = 10.0;
    double fastrand = 10.0; //Same as rand, we count their flops equally
    double floor = 1.0;
    double modulo = 15.0;
    double pow = 1.0; // 1 Mult  = pow2 (The only one we use)
    double exp = 13.0;
} throughputs_t;

static throughputs_t tp;

