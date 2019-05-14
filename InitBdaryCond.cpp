// #include <iostream>
// #include <stdlib.h>
// #include <math.h>
// #include <stdio.h>

#include "ConstPara.h"
// #include "matrix.h"
// #include "physvar.h"

matrix InitialSpace(matrix x)
{
    for (int i = 1; i <= N_In; i++)
    {
        x(1, i) = (i - 0.5) * dx;
    }
    return x;
}

physvar InitialCondition(matrix x)
{
    physvar U;
    double P;
    for (int i = 3; i <= N_In + 2; i++)
    {
        if (x(1, i-2) < 0.5 * L)
        {
            U.d(1, i) = 1.0;
            U.u(1, i) = 0.0;
            P = 1.0;
            U.E(1, i) = P / (gamma - 1.0) + 0.5 * U.d(1, i) * (U.u(1, i) * U.u(1, i));
        }
        else
        {
            U.d(1, i) = 0.125;
            U.u(1, i) = 0.0;
            P = 0.1;
            U.E(1, i) = P / (gamma - 1.0) + 0.5 * U.d(1, i) * (U.u(1, i) * U.u(1, i));
        }
    }
    U.u = dotprod(U.u, U.d);
    return U;
}


physvar BoundaryCondition(physvar U)
{// outflow
    for (int i = 1; i <= nghost; i++)
    {
        U.d(1, i) = U.d(1, nghost+1);
        U.u(1, i) = U.u(1, nghost+1);
        U.E(1, i) = U.E(1, nghost+1);

        U.d(1, N - i + 1) = U.d(1, N - nghost);
        U.u(1, N - i + 1) = U.u(1, N - nghost);
        U.E(1, N - i + 1) = U.E(1, N - nghost);
    }
    return U;
}