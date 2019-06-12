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

#include "Subfunc3D2std.cpp" // subfunction

#include <mpi.h>

using namespace std;


main(int argc, char *argv[])
{
    // initialize MPI
    int NRank, MyRank;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);

    MPI_Comm_size(MPI_COMM_WORLD, &NRank);

    // this test assumes only two ranks
    if (NRank != 2)
    {
        fprintf(stderr, "ERROR: NRank (%d) != 2\n", NRank);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    const int SendRank = 0;
    const int RecvRank = 1;
    // const int TargetRank = (MyRank + 1) % 2; // (0,1) --> (1,0) 
    const int NReq = 1;
    int Tag;
    MPI_Request Request;

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

    matrix fluxX_L(5);
    matrix fluxX_R(5);
    matrix dfluxX(5, N, NY, foo);
    matrix fluxX(5, N, NY, foo);
    matrix fluxY_L(5);
    matrix fluxY_R(5);
    matrix dfluxY(5, N, NY, foo);
    matrix fluxY(5, N, NY, foo);

    double dt;
    double *fluxPtr;
    double fluxtemp;
    double *UPtr;
    double Utemp;
    double SendBuf[N * NY * foo * 5] = {0.0};
    double RecvBuf[N * NY * foo * 5] = {0.0};

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

        if ( MyRank == RecvRank )
        {
            // data reconstruction
            DataReconstructionX_PLM(U, LX, RX);

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

            // compute fluxes
            // R[j-1] is the LEFT state at the j+1/2 inteface
            for (int m = nghost + 1; m <= N - nghost + 1; m++)
                for (int j = nghost + 1; j <= NY - nghost; j++)
                    for (int k = nghost + 1; k <= foo - nghost; k++)
                        fluxX.SetRow(m, j, k, RoeX(RX.GetRow(m - 1, j, k), LX.GetRow(m, j, k)));
        }
        else
        {
            // data reconstruction
            DataReconstructionY_PLM(U, LY, RY);

            // update the face-centered variables by 0.5*dt
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

            // compute fluxes
            // R[j-1] is the LEFT state at the j+1/2 inteface
            for (int m = nghost + 1; m <= N - nghost; m++)
                for (int j = nghost + 1; j <= NY - nghost + 1; j++)
                    for (int k = nghost + 1; k <= foo - nghost; k++)
                        fluxY.SetRow(m, j, k, RoeY(RY.GetRow(m, j - 1, k), LY.GetRow(m, j, k)));
        }

        Tag = 123;
        if (MyRank == SendRank)
        {
            for (int n = 0; n < 5; n++)
                for (int m = 0; m < N; m++)
                    for (int j = 0; j < NY; j++)
                        for (int k = 0; k < foo; k++)
                            SendBuf[k + foo * (j + NY * (m + N * n))] = fluxY(n + 1, m + 1, j + 1, k + 1);

            MPI_Isend(&SendBuf, N * NY * foo * 5, MPI_DOUBLE, RecvRank, Tag, MPI_COMM_WORLD, &Request);
        }
        else
        {
            MPI_Irecv(&RecvBuf, N * NY * foo * 5, MPI_DOUBLE, SendRank, Tag, MPI_COMM_WORLD, &Request);
        }
        // update time
        t = t + dt;

        // update the volume-averaged input variables by dt
        if (MyRank == RecvRank)
        {
            MPI_Wait(&Request, MPI_STATUSES_IGNORE);
            // printf("%f ", RecvBuf[5]);
            for (int n = 0; n < 5; n++)
                for (int m = 0; m < N; m++)
                    for (int j = 0; j < NY; j++)
                        for (int k = 0; k < foo; k++)
                            fluxY(n + 1, m + 1, j + 1, k + 1) = RecvBuf[k + foo * (j + NY * (m + N * n))];

            // fluxY.showY(10, 3);
            for (int m = nghost + 1; m <= N - nghost; m++)
                for (int j = nghost + 1; j <= NY - nghost; j++)
                    for (int k = nghost + 1; k <= foo - nghost; k++)
                    {
                        U.SetRow(m, j, k, U.GetRow(m, j, k) - dt / dx * (fluxX.GetRow(m + 1, j, k) - fluxX.GetRow(m, j, k)) - dt / dx * (fluxY.GetRow(m, j + 1, k) - fluxY.GetRow(m, j, k)));
                    }



            // write text
            if (myfile.is_open())
            {
                for (int m = nghost + 1; m <= N - nghost; m++)
                {
                    for (int j = nghost + 1; j <= N - nghost; j++)
                    {
                        // myfile << U(1, m, j, 3) << "\t";
                        // myfile << pow(U(2, m, j, 3) / U(1, m, j, 3), 2.0) + pow(U(3, m, j, 3) / U(1, m, j, 3), 2.0) << "\t";
                        myfile << ComputePressure(U(1, m, j, 3), U(2, m, j, 3), U(3, m, j, 3), U(4, m, j, 3), U(5, m, j, 3)) << "\t";
                    }
                    myfile << "\r\n";
                }
            }
            else
                cout << "Unable to open file";
        }

        Tag = 321;
        if (MyRank == RecvRank)
        {
            for (int n = 0; n < 5; n++)
                for (int m = 0; m < N; m++)
                    for (int j = 0; j < NY; j++)
                        for (int k = 0; k < foo; k++)
                            SendBuf[k + foo * (j + NY * (m + N * n))] = U(n + 1, m + 1, j + 1, k + 1);

            MPI_Isend(&SendBuf, N * NY * foo * 5, MPI_DOUBLE, SendRank, Tag, MPI_COMM_WORLD, &Request);
        }
        else
        {
            MPI_Irecv(&RecvBuf, N * NY * foo * 5, MPI_DOUBLE, RecvRank, Tag, MPI_COMM_WORLD, &Request);
        }        

        
        if (MyRank == SendRank)
        {
            MPI_Wait(&Request, MPI_STATUSES_IGNORE);

            for (int n = 0; n < 5; n++)
                for (int m = 0; m < N; m++)
                    for (int j = 0; j < NY; j++)
                        for (int k = 0; k < foo; k++)
                            U(n + 1, m + 1, j + 1, k + 1) = RecvBuf[k + foo * (j + NY * (m + N * n))];
        }
    }
    // terminate MPI
    MPI_Finalize();
    return 0;
}