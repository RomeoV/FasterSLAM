#pragma once
#include <cstdlib>
#include "flop_count.h"

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
    Instruction abs = AbsInstr(); // TODO
    Instruction sin = SinInstr(); //All values set to 1.0 for now, we change later
    Instruction atan2 = Atan2Instr();
    Instruction mul = MulInstr();
    Instruction add = AddInstr();
    Instruction div = DivInstr();
    Instruction sqrt = SqrtInstr();
    Instruction rsqrt = RsqrtInstr();
    Instruction negation = NegationInstr();
    Instruction cos = CosInstr();
    Instruction doublecomp = DoublecompInstr();
    Instruction rand = RandInstr(); //Guess
    Instruction fastrand = FastrandInstr(); //Guess
    Instruction floor = FloorInstr();
    Instruction modulo = ModuloInstr(); //Guess
    Instruction pow = PowInstr(); // 1 Mult  = pow2 (The only one we use)
    Instruction exp = ExpInstr();
} throughputs_t;

static throughputs_t tp;
