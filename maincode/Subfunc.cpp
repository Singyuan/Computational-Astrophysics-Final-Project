#ifndef __SUBFUNC_CPP__
#define __SUBFUNC_CPP__

#include "ConstPara.h"
#include "matrix.h"

// ##########################
// Define initial condition
// ##########################
matrix InitialCondition(double x)
{ // Sod shock tube
    matrix temp(1, 3);
    if (x < 0.5 * L)
    {
        temp(1, 1) = 1.0; // density
        temp(1, 2) = 0.0; // velocity x
        double P = 1.0;   // pressure
        temp(1, 3) = P / (mygamma - 1.0) + 0.5 * temp(1, 1) * (temp(1, 3) * temp(1, 3)); // energy density
    }
    else
    {
        temp(1, 1) = 0.125; // density
        temp(1, 2) = 0.0;   // velocity x
        double P = 0.1;     // pressure
        temp(1, 3) = P / (mygamma - 1.0) + 0.5 * temp(1, 1) * (temp(1, 3) * temp(1, 3)); // energy density
    }

    // conserved variables [0/1/2/3/4] <--> [density/momentum x/momentum y/momentum z/energy]
    temp(1, 2) = temp(1, 1) * temp(1, 2);
    return temp;
}

// ##########################
// Define boundary condition by setting ghost zones
// ##########################
void BoundaryCondition(matrix &U)
{ // outflow
    for (int i = 1; i <= nghost; i++)
    {
        U(i, 1) = U(nghost + 1, 1);
        U(i, 2) = U(nghost + 1, 2);
        U(i, 3) = U(nghost + 1, 3);

        U(N - i + 1, 1) = U(N - nghost, 1);
        U(N - i + 1, 2) = U(N - nghost, 2);
        U(N - i + 1, 3) = U(N - nghost, 3);
    }
}

// ##########################
// Compute pressure
// ##########################
double ComputePressure(double d, double px, double e)
{
    double P;
    P = (mygamma - 1.0) * (e - 0.5 * (pow(px, 2.0) / d));
    if (P <= 0)
    {
        
        printf("negative pressure in compute pressure !!\t%f\n", P);
        exit(1);
    }
        
    return P;
}

// ##########################
// Compute time-step by the CFL condition
// ##########################
double ComputeTimestep(matrix &U)
{
    double dt_cfl = 0.0;
    double max_info_speed = 0.0;
    double P, a, u;
    for (int i = 1; i <= N; i++)
    {
        P = ComputePressure(U(i, 1), U(i, 2), U(i, 3));
        a = pow((mygamma * P / U(i, 1)), 0.5);
        u = abs(U(i, 2) / U(i, 1));
        if (max_info_speed < u + a)
            max_info_speed = u + a;
    }
    dt_cfl = cfl * dx / max_info_speed;
    return dt_cfl;
}

// ##########################
// Compute limited slope
// ##########################
matrix ComputeLimitedSlope(matrix L, matrix C, matrix R)
{
    // compute the left and right slopes
    matrix slope_L(1, 3);
    matrix slope_R(1, 3);
    matrix slope_LR(1, 3);
    matrix slope_limited(1, 3);
    slope_L = C - L;
    slope_R = R - C;

    // apply the van-Leer limiter
    slope_LR = dotprod(slope_L, slope_R);
    for (int i = 1; i <= 3; i++)
    {
        if (slope_LR(1, i) > 0.0)
        {
            slope_limited(1, i) = 2.0 * slope_LR(1, i) / (slope_L(1, i) + slope_R(1, i));
        }
        else
        {
            slope_limited(1, i) = 0.0;
        }
    }
    return slope_limited;
}

// ##########################
// Convert conserved variables to primitive variables
// ##########################
matrix Conserved2Primitive(matrix U)
{
    matrix W(1, 3);
    W(1, 1) = U(1, 1);
    W(1, 2) = U(1, 2) / U(1, 1);
    W(1, 3) = ComputePressure(U(1, 1), U(1, 2), U(1, 3));
    return W;
}

// ##########################
// Convert primitive variables to conserved variables
// ##########################
matrix Primitive2Conserved(matrix W)
{
    matrix U(1, 3);
    U(1, 1) = W(1, 1);
    U(1, 2) = W(1, 1)*W(1, 2);
    U(1, 3) = W(1, 3) / (mygamma - 1.0) + 0.5 * W(1, 1)*pow(W(1, 2), 2.0);
    return U;
}

// ##########################
// Piecewise-linear data reconstruction (PLM)
// ##########################
void DataReconstruction_PLM(matrix &U, matrix &L, matrix &R)
{
    // allocate memory
    matrix W(N, 3);
    matrix slope(N, 3);

    // conserved variables-- > primitive variables
    for (int j = 1; j <= N; j++)
        W.SetRow(j, Conserved2Primitive(U.GetRow(j)));

    // compute the left and right states of each cell
    for (int j = 2; j <= N - 1; j++)
        slope.SetRow(j, ComputeLimitedSlope(W.GetRow(j - 1), W.GetRow(j), W.GetRow(j + 1)));

    matrix slope_limited(1, 3);
    for (int j = 2; j <= N - 1; j++)
    {
        //  get the face-centered variables
        slope_limited = slope.GetRow(j);
        L.SetRow(j, W.GetRow(j) - 0.5 * slope_limited);
        R.SetRow(j, W.GetRow(j) + 0.5 * slope_limited);

        // ensure face-centered variables lie between nearby volume-averaged (~cell-centered) values

        // primitive variables --> conserved variables
        L.SetRow(j, Primitive2Conserved(L.GetRow(j)));
        R.SetRow(j, Primitive2Conserved(R.GetRow(j)));
    }
}

// ##########################
// Piecewise-parabolic data reconstruction (PPM)
// ##########################
void DataReconstruction_PPM(matrix &U, matrix &L, matrix &R)
{
    // allocate memory
    matrix W(N, 3);
    matrix slope(N, 3);
    matrix acc(N, 3);

    // conserved variables-- > primitive variables
    for (int j = 1; j <= N; j++)
        W.SetRow(j, Conserved2Primitive(U.GetRow(j)));

    // compute the left and right states of each cell
    for (int j = 2; j <= N - 1; j++)
        slope.SetRow(j, ComputeLimitedSlope(W.GetRow(j - 1), W.GetRow(j), W.GetRow(j + 1)));

    // compute the left and right states of each cell
    for (int j = 3; j <= N - 2; j++)
        acc.SetRow(j, ComputeLimitedSlope(W.GetRow(j - 2), W.GetRow(j), W.GetRow(j + 2)));

    matrix slope_limited(1, 3);
    matrix acc_limited(1, 3);
    for (int j = 2; j <= N - 1; j++)
    {
        //  get the face-centered variables
        slope_limited = slope.GetRow(j);
        acc_limited = acc.GetRow(j);
        L.SetRow(j, W.GetRow(j) - 0.5 * slope_limited + 0.125 * acc_limited);
        R.SetRow(j, W.GetRow(j) + 0.5 * slope_limited + 0.125 * acc_limited);

        // ensure face-centered variables lie between nearby volume-averaged (~cell-centered) values

        // primitive variables --> conserved variables
        L.SetRow(j, Primitive2Conserved(L.GetRow(j)));
        R.SetRow(j, Primitive2Conserved(R.GetRow(j)));
    }
}

// ##########################
// Convert conserved variables to fluxes
// ##########################
matrix Conserved2Flux(matrix U)
{
    matrix flux(1, 3);
    double P = ComputePressure(U(1, 1), U(1, 2), U(1, 3));
    double v = U(1, 2) / U(1, 1);

    flux(1, 1) = U(1, 2);
    flux(1, 2) = v * U(1, 2) + P;
    flux(1, 3) = v * (U(1, 3) + P);
    return flux;
}

// ##########################
// Roe's Riemann solver
// ##########################
matrix Roe(matrix L, matrix R)
{
    // compute the enthalpy of the left and right states: H = (E+P)/rho
    matrix flux(1, 3);
    double P_L = ComputePressure(L(1, 1), L(1, 2), L(1, 3));
    double P_R = ComputePressure(R(1, 1), R(1, 2), R(1, 3));
    double H_L = (L(1, 3) + P_L) / L(1, 1);
    double H_R = (R(1, 3) + P_R) / R(1, 1);

    // compute Roe average values
    double rhoL_sqrt = pow(L(1, 1), 0.5);
    double rhoR_sqrt = pow(R(1, 1), 0.5);

    double u = (L(1, 2) / rhoL_sqrt + R(1, 2) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double H = (rhoL_sqrt * H_L + rhoR_sqrt * H_R) / (rhoL_sqrt + rhoR_sqrt);
    double V2 = u * u;
    double a = 0.0;

    // check negative pressure
    if (H - 0.5 * V2 > 0.0)
    {
        a = pow(((mygamma - 1.0) * (H - 0.5 * V2)), 0.5);
    }
    else
    {
        printf("negative pressure!\t%f\n", H - 0.5 * V2);
        exit(1);
    }

    // compute the amplitudes of different characteristic waves
    matrix dU(1, 3);
    matrix amp(1, 3);

    dU = R - L;
    amp(1, 2) = (mygamma - 1.0) / pow(a, 2.0) * (dU(1, 1) * (H - pow(u, 2.0)) + u * dU(1, 2) - dU(1, 3));
    amp(1, 1) = 0.5 / a * (dU(1, 1) * (u + a) - dU(1, 2) - a * amp(1, 2));
    amp(1, 3) = dU(1, 1) - amp(1, 1) - amp(1, 2);

    // compute the eigenvalues and right eigenvector matrix
    matrix EigVal(1, 3);
    matrix EigVec(3, 3);

    EigVal(1, 1) = abs(u - a);
    EigVal(1, 2) = abs(u);
    EigVal(1, 3) = abs(u + a);

    EigVec(1, 1) = 1.0;
    EigVec(1, 2) = u - a;
    EigVec(1, 3) = H - u * a;
    EigVec(2, 1) = 1.0;
    EigVec(2, 2) = u;
    EigVec(2, 3) = 0.5 * V2;
    EigVec(3, 1) = 1.0;
    EigVec(3, 2) = u + a;
    EigVec(3, 3) = H + u * a;

    // compute the fluxes of the left and right states
    matrix flux_L(1, 3);
    matrix flux_R(1, 3);

    flux_L = Conserved2Flux(L);
    flux_R = Conserved2Flux(R);

    // compute the Roe flux
    amp = dotprod(amp, EigVal);
    flux = 0.5 * (flux_L + flux_R) - 0.5 * amp * EigVec;
    return flux;
}


#endif //__SUBFUNC_CPP__