#ifndef __COMPUTEFLUX_CPP__
#define __COMPUTEFLUX_CPP__

#include "ConstPara.h"
#include "DataReconstruct.cpp"

// convert conserved variables to fluxes (overloading 1)
physvar Conserved2Flux(physvar U)
{
    physvar flux;
    matrix P(1, N);
    matrix v(1, N);
    P = ComputePressure(U.d, U.u, U.E);
    v = U.u / U.d;

    flux.d = U.u;
    flux.u = dotprod(v, U.u) + P;
    flux.E = dotprod(v, U.E + P);
    return flux;
}

// convert conserved variables to fluxes (overloading 2)
matrix Conserved2Flux(matrix U)
{
    matrix flux(3, 1);
    double P = ComputePressure(U(1, 1), U(2, 1), U(3, 1));
    double v = U(2, 1) / U(1, 1);

    flux(1, 1) = U(2, 1);
    flux(2, 1) = v*U(2, 1)+P;
    flux(3, 1) = v*(U(3, 1) + P);
    return flux;
}

// compute dflux the face-centered variables by 0.5*dt
matrix ComputedfluxSession(matrix FLd, matrix FRd, double dt)
{
    matrix dflux(1, N);
    for (int i = 2; i <= N -1; i++)
    {
        dflux(1, i) = 0.5 * dt / dx * (FRd(1, i) - FLd(1, i));
    }
    return dflux;
}

// Compute flux by Roe
matrix ComputeRoe(matrix L, matrix R)
{
    // compute the enthalpy of the left and right states: H = (E+P)/rho
    matrix flux(3, 1);
    double P_L = ComputePressure(L(1, 1), L(2, 1), L(3, 1));
    double P_R = ComputePressure(R(1, 1), R(2, 1), R(3, 1));
    double H_L = (L(3, 1) + P_L) / L(1, 1);
    double H_R = (R(3, 1) + P_R) / R(1, 1);

    // compute Roe average values
    double rhoL_sqrt = pow(L(1, 1), 0.5);
    double rhoR_sqrt = pow(R(1, 1), 0.5);

    double u = (L(2, 1) / rhoL_sqrt + R(2, 1) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double H = (rhoL_sqrt * H_L + rhoR_sqrt * H_R) / (rhoL_sqrt + rhoR_sqrt);
    double V2 = u * u;
    double a = 0;

    // check negative pressure
    if (H - 0.5 * V2 > 0.0)
    {
        a = pow(((gamma - 1.0) * (H - 0.5 * V2)), 0.5);
    }
    else
    {
        printf("negative pressure!\t%f\n", H - 0.5 * V2);
    }

    // compute the amplitudes of different characteristic waves
    matrix dU(3, 1);
    matrix amp(3, 1);

    dU = R - L;
    amp(2, 1) = (gamma - 1.0) / pow(a, 2.0) * (dU(1, 1) * (H - pow(u, 2.0)) + u * dU(2, 1) - dU(3, 1));
    amp(1, 1) = 0.5 / a * (dU(1, 1) * (u + a) - dU(2, 1) - a * amp(2, 1));
    amp(3, 1) = dU(1, 1) - amp(1, 1) - amp(2, 1);

    // compute the eigenvalues and right eigenvector matrix
    matrix EigVal(3, 1);
    matrix EigVec(3, 3);

    EigVal(1, 1) = abs(u - a);
    EigVal(2, 1) = abs(u);
    EigVal(3, 1) = abs(u + a);

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
    matrix flux_L(3, 1);
    matrix flux_R(3, 1);

    flux_L = Conserved2Flux(L);
    flux_R = Conserved2Flux(R);

    // compute the Roe flux
    amp = dotprod(amp, EigVal);
    EigVec = EigVec.T();
    flux = 0.5 * (flux_L + flux_R) - 0.5*EigVec*amp;
    return flux;
}

physvar Roe(physvar L, physvar R)
{
    matrix vectorL(3, 1); // temp matrix hold the vector on x
    matrix vectorR(3, 1);
    matrix tempflux(3, 1);
    physvar flux;

    for (int i = nghost + 1; i <= N - nghost; i++)
    {
        vectorL(1, 1) = L.d(1, i);
        vectorL(2, 1) = L.u(1, i);
        vectorL(3, 1) = L.E(1, i);

        vectorR(1, 1) = R.d(1, i - 1);
        vectorR(2, 1) = R.u(1, i - 1);
        vectorR(3, 1) = R.E(1, i - 1);

        tempflux = ComputeRoe(vectorR, vectorL);
        flux.d(1, i) = tempflux(1, 1);
        flux.u(1, i) = tempflux(2, 1);
        flux.E(1, i) = tempflux(3, 1);
    }
    return flux;
}


physvar ComputeNext(physvar flux, physvar U, double dt)
{
    for (int i = nghost + 1; i <= N - nghost;i++)
    {
        U.d(1, i) = U.d(1, i) - dt / dx * (flux.d(1, i)-flux.d(1, i-1));
        U.u(1, i) = U.u(1, i) - dt / dx * (flux.u(1, i) - flux.u(1, i - 1));
        U.E(1, i) = U.E(1, i) - dt / dx * (flux.E(1, i) - flux.E(1, i - 1));
    }
    return U;
}

#endif // __COMPUTEFLUX_CPP__