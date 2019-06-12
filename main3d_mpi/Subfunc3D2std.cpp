#ifndef __SUBFUNC_CPP__
#define __SUBFUNC_CPP__

#include "ConstPara.h"
#include "matrix3D.h"
// #include <math.h>
#include <cmath>
using namespace std;
// ##########################
// Define initial condition
// ##########################
matrix InitialCondition(double x, double y, double z)
{ // Sod shock tube
    matrix temp(5);
    if (sqrt(x * x + y * y + z * z) < 0.5 * L)
    {
        temp(1) = 1.0; // density
        temp(2) = 0.0; // velocity x
        temp(3) = 0.0; // velocity y
        temp(4) = 0.0; // velocity z
        double P = 1.0;   // pressure
        temp(5) = P / (mygamma - 1.0) + 0.5 * temp(1) * (pow(temp(2), 2.0) + pow(temp(3), 2.0) + pow(temp(4), 2.0)); // energy density
    }
    else
    {
        temp(1) = 0.125; // density
        temp(2) = 0.0;   // velocity x
        temp(3) = 0.0;   // velocity y
        temp(4) = 0.0;   // velocity z
        double P = 0.1;     // pressure
        temp(5) = P / (mygamma - 1.0) + 0.5 * temp(1) * (temp(3) * temp(3)); // energy density
    }

    // conserved variables [0/1/2/3/4] <--> [density/momentum x/momentum y/momentum z/energy]
    temp(2) = temp(1) * temp(2);
    temp(3) = temp(1) * temp(3);
    temp(4) = temp(1) * temp(4);
    return temp;
}

// ##########################
// Define boundary condition by setting ghost zones
// ##########################
void BoundaryCondition(matrix &U)
{   // outflow x
    for (int i = 1; i <= nghost; i++)
    {
        for (int j = 1; j <= NY_In; j++)
        {
            for (int k = 1; k <= foo_In; k++)
            {
                U(1, i, j + nghost, k + nghost) = U(1, nghost + 1, j + nghost, k + nghost);
                U(2, i, j + nghost, k + nghost) = U(2, nghost + 1, j + nghost, k + nghost);
                U(3, i, j + nghost, k + nghost) = U(3, nghost + 1, j + nghost, k + nghost);
                U(4, i, j + nghost, k + nghost) = U(4, nghost + 1, j + nghost, k + nghost);
                U(5, i, j + nghost, k + nghost) = U(5, nghost + 1, j + nghost, k + nghost);

                U(1, N - i + 1, j + nghost, k + nghost) = U(1, N - nghost, j + nghost, k + nghost);
                U(2, N - i + 1, j + nghost, k + nghost) = U(2, N - nghost, j + nghost, k + nghost);
                U(3, N - i + 1, j + nghost, k + nghost) = U(3, N - nghost, j + nghost, k + nghost);
                U(4, N - i + 1, j + nghost, k + nghost) = U(4, N - nghost, j + nghost, k + nghost);
                U(5, N - i + 1, j + nghost, k + nghost) = U(5, N - nghost, j + nghost, k + nghost);
            }
        }   
    }

    // outflow y
    for (int i = 1; i <= N_In; i++)
    {
        for (int j = 1; j <= nghost; j++)
        {
            for (int k = 1; k <= foo_In; k++)
            {
                U(1, i + nghost, j, k + nghost) = U(1, i + nghost, nghost + 1, k + nghost);
                U(2, i + nghost, j, k + nghost) = U(2, i + nghost, nghost + 1, k + nghost);
                U(3, i + nghost, j, k + nghost) = U(3, i + nghost, nghost + 1, k + nghost);
                U(4, i + nghost, j, k + nghost) = U(4, i + nghost, nghost + 1, k + nghost);
                U(5, i + nghost, j, k + nghost) = U(5, i + nghost, nghost + 1, k + nghost);

                U(1, i + nghost, N - j + 1, k + nghost) = U(1, i + nghost, N - nghost, k + nghost);
                U(2, i + nghost, N - j + 1, k + nghost) = U(2, i + nghost, N - nghost, k + nghost);
                U(3, i + nghost, N - j + 1, k + nghost) = U(3, i + nghost, N - nghost, k + nghost);
                U(4, i + nghost, N - j + 1, k + nghost) = U(4, i + nghost, N - nghost, k + nghost);
                U(5, i + nghost, N - j + 1, k + nghost) = U(5, i + nghost, N - nghost, k + nghost);
            }
        }
    }

    // // outflow z
    // for (int i = 1; i <= N_In; i++)
    // {
    //     for (int j = 1; j <= NY_In; j++)
    //     {
    //         for (int k = 1; k <= nghost; k++)
    //         {
    //             U(1, i + nghost, j + nghost, k) = U(1, i + nghost, j + nghost, nghost + 1);
    //             U(2, i + nghost, j + nghost, k) = U(2, i + nghost, j + nghost, nghost + 1);
    //             U(3, i + nghost, j + nghost, k) = U(3, i + nghost, j + nghost, nghost + 1);
    //             U(4, i + nghost, j + nghost, k) = U(4, i + nghost, j + nghost, nghost + 1);
    //             U(5, i + nghost, j + nghost, k) = U(5, i + nghost, j + nghost, nghost + 1);

    //             U(1, i + nghost, j + nghost, N - k + 1) = U(1, i + nghost, j + nghost, N - nghost);
    //             U(2, i + nghost, j + nghost, N - k + 1) = U(2, i + nghost, j + nghost, N - nghost);
    //             U(3, i + nghost, j + nghost, N - k + 1) = U(3, i + nghost, j + nghost, N - nghost);
    //             U(4, i + nghost, j + nghost, N - k + 1) = U(4, i + nghost, j + nghost, N - nghost);
    //             U(5, i + nghost, j + nghost, N - k + 1) = U(5, i + nghost, j + nghost, N - nghost);
    //         }
    //     }
    // }
}

// ##########################
// Compute pressure
// ##########################
double ComputePressure(double d, double px, double py, double pz, double e)
{
    double P;
    P = (mygamma - 1.0) * (e - 0.5 * (pow(px, 2.0) + pow(py, 2.0) + pow(pz, 2.0) / d));
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
    double P, a, u, v, w;
    for (int i = 1; i <= N; i++)
    {
        for (int j = 1; j <= NY; j++)
        {
            for (int k = 1; k <= foo; k++)
            {
                P = ComputePressure(U(1, i, j, k), U(2, i, j, k), U(3, i, j, k), U(4, i, j, k), U(5, i, j, k));
                a = pow((mygamma * P / U(1, i, j, k)), 0.5);
                u = std::abs(U(2, i, j, k) / U(1, i, j, k));
                v = std::abs(U(3, i, j, k) / U(1, i, j, k));
                w = std::abs(U(4, i, j, k) / U(1, i, j, k));
                if (max_info_speed < u + a)
                    max_info_speed = u + a;
                if (max_info_speed < v + a)
                    max_info_speed = v + a;
                // if (max_info_speed < w + a)
                //     max_info_speed = w + a;
            }
        }
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
    matrix slope_L(5);
    matrix slope_R(5);
    matrix slope_LR(5);
    matrix slope_limited(5);
    slope_L = C - L;
    slope_R = R - C;

    // apply the van-Leer limiter
    slope_LR = dotprod(slope_L, slope_R);
    for (int i = 1; i <= 5; i++)
    {
        if (slope_LR(i) > 0.0)
        {
            slope_limited(i) = 2.0 * slope_LR(i) / (slope_L(i) + slope_R(i));
        }
        else
        {
            slope_limited(i) = 0.0;
        }
    }
    return slope_limited;
}

// ##########################
// Convert conserved variables to primitive variables
// ##########################
matrix Conserved2Primitive(matrix U)
{
    matrix W(5);
    W(1) = U(1);
    W(2) = U(2) / U(1);
    W(3) = U(3) / U(1);
    W(4) = U(4) / U(1);
    W(5) = ComputePressure(U(1), U(2), U(3), U(4), U(5));
    return W;
}

// ##########################
// Convert primitive variables to conserved variables
// ##########################
matrix Primitive2Conserved(matrix W)
{
    matrix U(5);
    U(1) = W(1);
    U(2) = W(1) * W(2);
    U(3) = W(1) * W(3);
    U(4) = W(1) * W(4);
    U(5) = W(5) / (mygamma - 1.0) + 0.5 * W(1) * (pow(W(2), 2.0) + pow(W(3), 2.0) + pow(W(4), 2.0));
    return U;
}

// ##########################
// Piecewise-linear data reconstruction (PLM)
// ##########################
void DataReconstructionX_PLM(matrix &U, matrix &L, matrix &R)
{
    // allocate memory
    matrix W(5, N, NY, foo);
    matrix slope(5, N, NY, foo);

    // conserved variables-- > primitive variables
    for (int i = 1; i <= N; i++)
        for (int j = nghost + 1; j <= NY - nghost; j++)
            for (int k = nghost + 1; k <= foo - nghost; k++)
                W.SetRow(i, j, k, Conserved2Primitive(U.GetRow(i, j, k)));

    // compute the left and right states of each cell
    for (int i = 2; i <= N - 1; i++)
        for (int j = nghost + 1; j <= NY - nghost; j++)
            for (int k = nghost + 1; k <= foo - nghost; k++)
                slope.SetRow(i, j, k, ComputeLimitedSlope(W.GetRow(i - 1, j, k), W.GetRow(i, j, k), W.GetRow(i + 1, j, k)));

    matrix slope_limited(5);
    for (int i = 2; i <= N - 1; i++)
    {
        for (int j = nghost + 1; j <= NY - nghost; j++)
        {
            for (int k = nghost + 1; k <= foo - nghost; k++)
            {
                //  get the face-centered variables
                slope_limited = slope.GetRow(i, j, k);
                L.SetRow(i, j, k, W.GetRow(i, j, k) - 0.5 * slope_limited);
                R.SetRow(i, j, k, W.GetRow(i, j, k) + 0.5 * slope_limited);

                // ensure face-centered variables lie between nearby volume-averaged (~cell-centered) values
                for (int n = 1; n <= 5; n++)
                {
                    L(n, i, j, k) = max(L(n, i, j, k), min(W(n, i - 1, j, k), W(n, i, j, k)));
                    L(n, i, j, k) = min(L(n, i, j, k), max(W(n, i - 1, j, k), W(n, i, j, k)));
                    R(n, i, j, k) = 2.0 * W(n, i, j, k) - L(n, i, j, k);

                    R(n, i, j, k) = max(R(n, i, j, k), min(W(n, i + 1, j, k), W(n, i, j, k)));
                    R(n, i, j, k) = min(R(n, i, j, k), max(W(n, i + 1, j, k), W(n, i, j, k)));
                    L(n, i, j, k) = 2.0 * W(n, i, j, k) - R(n, i, j, k);
                }

                // primitive variables --> conserved variables
                L.SetRow(i, j, k, Primitive2Conserved(L.GetRow(i, j, k)));
                R.SetRow(i, j, k, Primitive2Conserved(R.GetRow(i, j, k)));
            }
        }
    }
}

void DataReconstructionY_PLM(matrix &U, matrix &L, matrix &R)
{
    // allocate memory
    matrix W(5, N, NY, foo);
    matrix slope(5, N, NY, foo);

    // conserved variables-- > primitive variables
    for (int i = nghost + 1; i <= N - nghost; i++)
        for (int j = 1; j <= NY; j++)
            for (int k = nghost + 1; k <= foo - nghost; k++)
                W.SetRow(i, j, k, Conserved2Primitive(U.GetRow(i, j, k)));

    // compute the left and right states of each cell
    for (int i = nghost + 1; i <= N - nghost; i++)
        for (int j = 2; j <= NY - 1; j++)
            for (int k = nghost + 1; k <= foo - nghost; k++)
                slope.SetRow(i, j, k, ComputeLimitedSlope(W.GetRow(i, j - 1, k), W.GetRow(i, j, k), W.GetRow(i, j + 1, k)));

    matrix slope_limited(5);
    for (int i = nghost + 1; i <= N - nghost; i++)
    {
        for (int j = 2; j <= NY - 1; j++)
        {
            for (int k = nghost + 1; k <= foo - nghost; k++)
            {
                //  get the face-centered variables
                slope_limited = slope.GetRow(i, j, k);
                L.SetRow(i, j, k, W.GetRow(i, j, k) - 0.5 * slope_limited);
                R.SetRow(i, j, k, W.GetRow(i, j, k) + 0.5 * slope_limited);

                // ensure face-centered variables lie between nearby volume-averaged (~cell-centered) values
                for (int n = 1; n <= 5; n++)
                {
                    L(n, i, j, k) = max(L(n, i, j, k), min(W(n, i, j - 1, k), W(n, i, j, k)));
                    L(n, i, j, k) = min(L(n, i, j, k), max(W(n, i, j - 1, k), W(n, i, j, k)));
                    R(n, i, j, k) = 2.0 * W(n, i, j, k) - L(n, i, j, k);

                    R(n, i, j, k) = max(R(n, i, j, k), min(W(n, i, j + 1, k), W(n, i, j, k)));
                    R(n, i, j, k) = min(R(n, i, j, k), max(W(n, i, j + 1, k), W(n, i, j, k)));
                    L(n, i, j, k) = 2.0 * W(n, i, j, k) - R(n, i, j, k);
                }
                // primitive variables --> conserved variables
                L.SetRow(i, j, k, Primitive2Conserved(L.GetRow(i, j, k)));
                R.SetRow(i, j, k, Primitive2Conserved(R.GetRow(i, j, k)));
            }
        }
    }
}

void DataReconstructionZ_PLM(matrix &U, matrix &L, matrix &R)
{
    // allocate memory
    matrix W(5, N, NY, foo);
    matrix slope(5, N, NY, foo);

    // conserved variables-- > primitive variables
    for (int i = nghost + 1; i <= N - nghost; i++)
        for (int j = nghost + 1; j <= NY - nghost; j++)
            for (int k = 1; k <= foo; k++)
                W.SetRow(i, j, k, Conserved2Primitive(U.GetRow(i, j, k)));

    // compute the left and right states of each cell
    for (int i = nghost + 1; i <= N - nghost; i++)
        for (int j = nghost + 1; j <= NY - nghost; j++)
            for (int k = 2; k <= foo - 1; k++)
                slope.SetRow(i, j, k, ComputeLimitedSlope(W.GetRow(i, j, k - 1), W.GetRow(i, j, k), W.GetRow(i, j, k + 1)));

    matrix slope_limited(5);
    for (int i = nghost + 1; i <= N - nghost; i++)
    {
        for (int j = nghost + 1; j <= NY - nghost; j++)
        {
            for (int k = 2; k <= foo - 1; k++)
            {
                //  get the face-centered variables
                slope_limited = slope.GetRow(i, j, k);
                L.SetRow(i, j, k, W.GetRow(i, j, k) - 0.5 * slope_limited);
                R.SetRow(i, j, k, W.GetRow(i, j, k) + 0.5 * slope_limited);

                // ensure face-centered variables lie between nearby volume-averaged (~cell-centered) values

                // primitive variables --> conserved variables
                L.SetRow(i, j, k, Primitive2Conserved(L.GetRow(i, j, k)));
                R.SetRow(i, j, k, Primitive2Conserved(R.GetRow(i, j, k)));
            }
        }
    }
    // L.display();
}

// // ##########################
// // Piecewise-parabolic data reconstruction (PPM)
// // ##########################
// void DataReconstruction_PPM(matrix &U, matrix &L, matrix &R)
// {
//     // allocate memory
//     matrix W(N, 3);
//     matrix slope(N, 3);
//     matrix acc(N, 3);

//     // conserved variables-- > primitive variables
//     for (int j = 1; j <= N; j++)
//         W.SetRow(j, Conserved2Primitive(U.GetRow(j)));

//     // compute the left and right states of each cell
//     for (int j = 2; j <= N - 1; j++)
//         slope.SetRow(j, ComputeLimitedSlope(W.GetRow(j - 1), W.GetRow(j), W.GetRow(j + 1)));

//     // compute the left and right states of each cell
//     for (int j = 3; j <= N - 2; j++)
//         acc.SetRow(j, ComputeLimitedSlope(slope.GetRow(j - 1), slope.GetRow(j), slope.GetRow(j + 1)));

//     matrix slope_limited(1, 3);
//     matrix acc_limited(1, 3);
//     for (int j = 2; j <= N - 1; j++)
//     {
//         //  get the face-centered variables
//         slope_limited = slope.GetRow(j);
//         acc_limited = acc.GetRow(j);
//         L.SetRow(j, W.GetRow(j) - 0.5 * slope_limited + 0.125 * acc_limited);
//         R.SetRow(j, W.GetRow(j) + 0.5 * slope_limited + 0.125 * acc_limited);

//         // ensure face-centered variables lie between nearby volume-averaged (~cell-centered) values

//         // primitive variables --> conserved variables
//         L.SetRow(j, Primitive2Conserved(L.GetRow(j)));
//         R.SetRow(j, Primitive2Conserved(R.GetRow(j)));
//     }
// }

// ##########################
// Convert conserved variables to fluxes
// ##########################
matrix Conserved2FluxX(matrix U)
{
    matrix flux(5);
    double P = ComputePressure(U(1), U(2), U(3), U(4), U(5));
    double v = U(2) / U(1);
    flux(1) = U(2);
    flux(2) = v * U(2) + P;
    flux(3) = v * U(3);
    flux(4) = v * U(4);
    flux(5) = v * (U(5) + P);
    return flux;
}

matrix Conserved2FluxY(matrix U)
{
    matrix flux(5);
    double P = ComputePressure(U(1), U(2), U(3), U(4), U(5));
    double v = U(3) / U(1);
    flux(1) = U(3);
    flux(2) = v * U(2);
    flux(3) = v * U(3) + P;
    flux(4) = v * U(4);
    flux(5) = v * (U(5) + P);
    return flux;
}

matrix Conserved2FluxZ(matrix U)
{
    matrix flux(5);
    double P = ComputePressure(U(1), U(2), U(3), U(4), U(5));
    double v = U(4) / U(1);
    flux(1) = U(4);
    flux(2) = v * U(3);
    flux(3) = v * U(3);
    flux(4) = v * U(4) + P;
    flux(5) = v * (U(5) + P);
    return flux;
}

matrix Conserved2FluxAddressX(matrix &U)
{
    matrix flux(5);
    double P = ComputePressure(U(1), U(2), U(3), U(4), U(5));
    double v = U(2) / U(1);
    flux(1) = U(2);
    flux(2) = v * U(2) + P;
    flux(3) = v * U(3);
    flux(4) = v * U(4);
    flux(5) = v * (U(5) + P);
    return flux;
}

matrix Conserved2FluxAddressY(matrix &U)
{
    matrix flux(5);
    double P = ComputePressure(U(1), U(2), U(3), U(4), U(5));
    double v = U(3) / U(1);
    flux(1) = U(2);
    flux(2) = v * U(2);
    flux(3) = v * U(3) + P;
    flux(4) = v * U(4);
    flux(5) = v * (U(5) + P);
    return flux;
}

matrix Conserved2FluxAddressZ(matrix &U)
{
    matrix flux(5);
    double P = ComputePressure(U(1), U(2), U(3), U(4), U(5));
    double v = U(4) / U(1);
    flux(1) = U(2);
    flux(2) = v * U(2);
    flux(3) = v * U(3);
    flux(4) = v * U(4) + P;
    flux(5) = v * (U(5) + P);
    return flux;
}

// ##########################
// Roe's Riemann solver
// ##########################
matrix RoeX(matrix L, matrix R)
{
    // compute the enthalpy of the left and right states: H = (E+P)/rho
    matrix flux(5);
    double P_L = ComputePressure(L(1), L(2), L(3), L(4), L(5));
    double P_R = ComputePressure(R(1), R(2), R(3), R(4), R(5));
    double H_L = (L(5) + P_L) / L(1);
    double H_R = (R(5) + P_R) / R(1);
    
    // compute Roe average values
    double rhoL_sqrt = pow(L(1), 0.5);
    double rhoR_sqrt = pow(R(1), 0.5);

    double u = (L(2) / rhoL_sqrt + R(2) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double v = (L(3) / rhoL_sqrt + R(3) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double w = (L(4) / rhoL_sqrt + R(4) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double H = (rhoL_sqrt * H_L + rhoR_sqrt * H_R) / (rhoL_sqrt + rhoR_sqrt);
    double V2 = u * u + v * v + w * w;
    double a = 0.0;

    // check negative pressure
    if (H - 0.5 * V2 > 0.0)
    {
        a = pow(((mygamma - 1.0) * (H - 0.5 * V2)), 0.5);
    }
    else
    {
        printf("negative pressure in Roe !\t%f\n", H - 0.5 * V2);
        exit(1);
    }

    // compute the amplitudes of different characteristic waves
    matrix dU(5);
    matrix amp(5);

    dU = R - L;
    amp(3) = dU(3) - v * dU(1);
    amp(4) = dU(4) - w * dU(1);
    amp(2) = (mygamma - 1.0) / pow(a, 2.0) * (dU(1) * (H - pow(u, 2.0)) + u * dU(2) - dU(5) + v * amp(3) + w * amp(4));
    amp(1) = 0.5 / a * (dU(1, 1) * (u + a) - dU(2) - a * amp(2));
    amp(5) = dU(1) - amp(1) - amp(2);

    // compute the eigenvalues and right eigenvector matrix
    matrix EigVal(5);
    matrix EigVec(5, 5);

    EigVal(1) = std::abs(u - a);
    EigVal(2) = std::abs(u);
    EigVal(3) = std::abs(u);
    EigVal(4) = std::abs(u);
    EigVal(5) = std::abs(u + a);

    EigVec(1, 1) = 1.0;
    EigVec(2, 1) = u - a;
    EigVec(3, 1) = v;
    EigVec(4, 1) = w;
    EigVec(5, 1) = H - u * a;
    EigVec(1, 2) = 1.0;
    EigVec(2, 2) = u;
    EigVec(3, 2) = v;
    EigVec(4, 2) = w;
    EigVec(5, 2) = 0.5 * V2;
    EigVec(3, 3) = 1.0;
    EigVec(5, 3) = v;
    EigVec(4, 4) = 1.0;
    EigVec(5, 4) = w;
    EigVec(1, 5) = 1.0;
    EigVec(2, 5) = u + a;
    EigVec(3, 5) = v;
    EigVec(4, 5) = w;
    EigVec(5, 5) = H + u * a;

    // compute the fluxes of the left and right states
    matrix flux_L(5);
    matrix flux_R(5);
    flux_L = Conserved2FluxAddressX(L);
    flux_R = Conserved2FluxAddressX(R);

    // // compute the Roe flux
    amp = dotprod(amp, EigVal);
    flux = 0.5 * (flux_L + flux_R) - 0.5 * EigVec * amp;
    return flux;
}

matrix RoeY(matrix L, matrix R)
{
    // compute the enthalpy of the left and right states: H = (E+P)/rho
    matrix flux(5);
    double P_L = ComputePressure(L(1), L(2), L(3), L(4), L(5));
    double P_R = ComputePressure(R(1), R(2), R(3), R(4), R(5));
    double H_L = (L(5) + P_L) / L(1);
    double H_R = (R(5) + P_R) / R(1);

    // compute Roe average values
    double rhoL_sqrt = pow(L(1), 0.5);
    double rhoR_sqrt = pow(R(1), 0.5);

    double u = (L(2) / rhoL_sqrt + R(2) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double v = (L(3) / rhoL_sqrt + R(3) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double w = (L(4) / rhoL_sqrt + R(4) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double H = (rhoL_sqrt * H_L + rhoR_sqrt * H_R) / (rhoL_sqrt + rhoR_sqrt);
    double V2 = u * u + v * v + w * w;
    double a = 0.0;
    
    // check negative pressure
    if (H - 0.5 * V2 > 0.0)
    {
        a = pow(((mygamma - 1.0) * (H - 0.5 * V2)), 0.5);
    }
    else
    {
        printf("negative pressure in Roe !\t%f\n", H - 0.5 * V2);
        exit(1);
    }

    // compute the amplitudes of different characteristic waves
    matrix dU(5);
    matrix amp(5);

    dU = R - L;
    amp(2) = dU(2) - u * dU(1);
    amp(4) = dU(4) - w * dU(1);
    amp(3) = (mygamma - 1.0) / pow(a, 2.0) * (dU(1) * (H - pow(v, 2.0)) + v * dU(3) - dU(5) + u * amp(2) + w * amp(4));
    amp(1) = 0.5 / a * (dU(1, 1) * (v + a) - dU(3) - a * amp(3));
    amp(5) = dU(1) - amp(1) - amp(3);

    // compute the eigenvalues and right eigenvector matrix
    matrix EigVal(5);
    matrix EigVec(5, 5);

    EigVal(1) = std::abs(v - a);
    EigVal(2) = std::abs(v);
    EigVal(3) = std::abs(v);
    EigVal(4) = std::abs(v);
    EigVal(5) = std::abs(v + a);

    EigVec(1, 1) = 1.0;
    EigVec(2, 1) = u;
    EigVec(3, 1) = v - a;
    EigVec(4, 1) = w;
    EigVec(5, 1) = H - v * a;
    EigVec(2, 2) = 1.0;
    EigVec(5, 2) = u;
    EigVec(1, 3) = 1.0;
    EigVec(2, 3) = u;
    EigVec(3, 3) = v;
    EigVec(4, 3) = w;
    EigVec(5, 3) = 0.5 * V2;
    EigVec(4, 4) = 1.0;
    EigVec(5, 4) = w;
    EigVec(1, 5) = 1.0;
    EigVec(2, 5) = u;
    EigVec(3, 5) = v + a;
    EigVec(4, 5) = w;
    EigVec(5, 5) = H + v * a;

    // compute the fluxes of the left and right states
    matrix flux_L(5);
    matrix flux_R(5);
    flux_L = Conserved2FluxAddressY(L);
    flux_R = Conserved2FluxAddressY(R);

    // // compute the Roe flux
    amp = dotprod(amp, EigVal);
    flux = 0.5 * (flux_L + flux_R) - 0.5 * EigVec * amp;
    return flux;
}



#endif //__SUBFUNC_CPP__