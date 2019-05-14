#ifndef __SUBFUNC_CPP__
#define __SUBFUNC_CPP__

#include "ConstPara.h"

// compute pressure (overloading 1)
matrix ComputePressure(matrix d, matrix px, matrix e)
{
    matrix P(1, N);
    P = (gamma - 1.0) * (e - 0.5 * (dotpow(px, 2.0) / d));
    return P;
}

// compute pressure (overloading 2)
double ComputePressure(double d, double px, double e)
{
    double P;
    P = (gamma - 1.0) * (e - 0.5 * (pow(px, 2.0) / d));
    return P;
}

// convert conserved variables to primitive variables
physvar Conserved2Primitive(physvar U)
{
    physvar W;
    W.d = U.d;
    W.u = U.u / U.d;
    W.E = ComputePressure(U.d, U.u, U.E);
    return W;
}

// convert primitive variables to conserved variables
physvar Primitive2Conserved(physvar W)
{
    physvar U;
    U.d = W.d;
    U.u = dotprod(W.u, W.d);
    U.E = W.E / (gamma - 1.0) + 0.5 * dotprod(W.d, dotpow(W.u, 2));
    return U;
}

#endif //__SUBFUNC_CPP__