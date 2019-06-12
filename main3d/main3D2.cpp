#include <iostream>
// #include <stdio.h>
// #include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <string.h>
// #include <time.h>
#include <fstream> // text
#include <omp.h> // use open mp

#include "ConstPara.h" // constant parameter
#include "matrix3D.h" // omp version

#include "Subfunc3D2.cpp" // subfunction

using namespace std;


main(int argc, char *argv[])
{
    // set open mp
    // omp_set_num_threads(NThread);

    // set initial condition
    double t = 0.0;     // time axis
    matrix x(N_In);     // space axis
    matrix U(5, N, NY, foo);     // physics variable
    for (int i = 1; i <= N_In; i++)
    {
        x(i) = (i - 0.5) * dx; // cell-centered coordinates
    }

    for (int i = 1; i <= N_In; i++)
    {
        for (int j = 1; j <= NY_In; j++)
        {
            for (int k = 1; k <= foo_In; k++)
            {
                U.SetRow(i + nghost, j + nghost, k + nghost, InitialCondition(x(i), x(j), 0));
            }
        }
    }

    // set flux variable hold the memory
    matrix LX(5, N, NY, foo);
    matrix RX(5, N, NY, foo);
    matrix LY(5, N, NY, foo);
    matrix RY(5, N, NY, foo);
    matrix LZ(5, N, NY, foo);
    matrix RZ(5, N, NY, foo);

    matrix fluxX_L(5);
    matrix fluxX_R(5);
    matrix dfluxX(5, N, NY, foo);
    matrix fluxX(5, N, NY, foo);
    matrix fluxY_L(5);
    matrix fluxY_R(5);
    matrix dfluxY(5, N, NY, foo);
    matrix fluxY(5, N, NY, foo);
    // matrix fluxZ_L(5);
    // matrix fluxZ_R(5);
    // matrix dfluxZ(5, N, NY, foo);
    // matrix fluxZ(5, N, NY, foo);
    double dt;

    // open text
    ofstream myfile("anidata.txt");
    for (int i = 1; i <= N_In; i++)
    {
        myfile << x(i) << "\t";
    }
    myfile << "\r\n";

    // main loop
    for (int i = 0; i < largenumber; i++)
    {
        // set the boundary conditions
        BoundaryCondition(U);

        // estimate time-step from the CFL condition
        dt = ComputeTimestep(U);
        // dt = 0.002;
        printf("t = %3.4f --> %3.4f\n", t, t + dt);

        // data reconstruction
        DataReconstructionX_PLM(U, LX, RX);
        DataReconstructionY_PLM(U, LY, RY);
        // DataReconstructionZ_PLM(U, LZ, RZ);

        // update the face-centered variables by 0.5*dt
        for (int m = 2; m <= N - 1; m++)
        {
            for (int j = nghost + 1; j <= NY - nghost; j++)
            {
                for (int k = nghost + 1; k <= foo - nghost; k++)
                {
                    fluxX_L = Conserved2FluxX(LX.GetRow(m, j, k));
                    fluxX_R = Conserved2FluxX(RX.GetRow(m, j, k));
                    dfluxX.SetRow(m, j, k, 0.5 * dt / dx * (fluxX_R - fluxX_L));
                    LX.SetRow(m, j, k, LX.GetRow(m, j, k) - dfluxX.GetRow(m, j, k));
                    RX.SetRow(m, j, k, RX.GetRow(m, j, k) - dfluxX.GetRow(m, j, k));
                }
            }
        }
        

        for (int m = nghost + 1; m <= N - nghost; m++)
        {
            for (int j = 2; j <= NY - 1; j++)
            {
                for (int k = nghost + 1; k <= foo - nghost; k++)
                {
                    fluxY_L = Conserved2FluxY(LY.GetRow(m, j, k));
                    fluxY_R = Conserved2FluxY(RY.GetRow(m, j, k));
                    dfluxY.SetRow(m, j, k, 0.5 * dt / dx * (fluxY_R - fluxY_L));
                    LY.SetRow(m, j, k, LY.GetRow(m, j, k) - dfluxY.GetRow(m, j, k));
                    RY.SetRow(m, j, k, RY.GetRow(m, j, k) - dfluxY.GetRow(m, j, k));
                }
            }
        }
        
        // for (int m = nghost + 1; m <= N - nghost; m++)
        // {
        //     for (int j = nghost + 1; j <= NY - nghost; j++)
        //     {
        //         for (int k = 2; k <= foo - 1; k++)
        //         {
        //             fluxZ_L = Conserved2FluxZ(RZ.GetRow(m, j, k));
        //             fluxZ_R = Conserved2FluxZ(RZ.GetRow(m, j, k));
        //             dfluxZ.SetRow(m, j, k, 0.5 * dt / dx * (fluxZ_R - fluxZ_L));
        //             LZ.SetRow(m, j, k, LZ.GetRow(m, j, k) - dfluxZ.GetRow(m, j, k));
        //             RZ.SetRow(m, j, k, RZ.GetRow(m, j, k) - dfluxZ.GetRow(m, j, k));
        //         }
        //     }
        // }

        // compute fluxes
        // R[j-1] is the LEFT state at the j+1/2 inteface
        for (int m = nghost + 1; m <= N - nghost + 1; m++)
            for (int j = nghost + 1; j <= NY - nghost; j++)
                for (int k = nghost + 1; k <= foo - nghost; k++)
                    fluxX.SetRow(m, j, k, RoeX(RX.GetRow(m - 1, j, k), LX.GetRow(m, j, k)));

        for (int m = nghost + 1; m <= N - nghost; m++)
            for (int j = nghost + 1; j <= NY - nghost + 1; j++)
                for (int k = nghost + 1; k <= foo - nghost; k++)
                    fluxY.SetRow(m, j, k, RoeY(RY.GetRow(m, j - 1, k), LY.GetRow(m, j, k)));

        // for (int m = nghost + 1; m <= N - nghost; m++)
        //     for (int j = nghost + 1; j <= NY - nghost; j++)
        //         for (int k = nghost + 1; k <= foo - nghost + 1; k++)
        //             fluxZ.SetRow(m, j, k, RoeZ(RZ.GetRow(m, j, k - 1), LZ.GetRow(m, j, k)));

        // update the volume-averaged input variables by dt
        for (int m = nghost + 1; m <= N - nghost; m++)
            for (int j = nghost + 1; j <= NY - nghost; j++)
                for (int k = nghost + 1; k <= foo - nghost; k++)
                {
                    // U.SetRow(m, j, k, U.GetRow(m, j, k) - dt / dx * (fluxX.GetRow(m + 1, j, k) - fluxX.GetRow(m, j, k)) - dt / dx * (fluxY.GetRow(m, j + 1, k) - fluxY.GetRow(m, j, k)) - dt / dx *(fluxZ.GetRow(m, j, k + 1) - fluxZ.GetRow(m, j, k)));
                    U.SetRow(m, j, k, U.GetRow(m, j, k) - dt / dx * (fluxX.GetRow(m + 1, j, k) - fluxX.GetRow(m, j, k)) - dt / dx * (fluxY.GetRow(m, j + 1, k) - fluxY.GetRow(m, j, k)));
                }

        // update time
        t = t + dt;

        // write text
        if (myfile.is_open())
        {
            for (int m = nghost + 1; m <= N - nghost; m++)
            {
                for (int j = nghost + 1; j <= N - nghost; j++)
                {
                    // myfile << U(1, m, j, 3) << "\t";
                    // myfile << sqrt(pow(U(2, m, j, 3) / U(1, m, j, 3), 2.0) + pow(U(3, m, j, 3) / U(1, m, j, 3), 2.0)) << "\t";
                    myfile << ComputePressure(U(1, m, j, 3), U(2, m, j, 3), U(3, m, j, 3), U(4, m, j, 3), U(5, m, j, 3)) << "\t";
                }
                myfile << "\r\n";
            }
        }
        else
            cout << "Unable to open file";
    }
    return 0;
}