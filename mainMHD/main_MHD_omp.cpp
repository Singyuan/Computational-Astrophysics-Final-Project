#include <iostream>
// #include <stdio.h>
// #include <stdlib.h>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <string.h>
#include <fstream> // text
#include <omp.h> // use open mp
#include <algorithm>
#include "ConstPara.h" // constant parameter
#include "matrix.h" // omp version

#include "Subfunc_MHD.cpp" // subfunction

using namespace std;




main(int argc, char *argv[])
{
    // set open mp
    omp_set_num_threads(NThread);

    // set initial condition
    double t = 0.0;     // time axis
    matrix x(1, N_In);  // space axis
    matrix U(N, 8);     // physics variable
    for (int i = 1; i <= N_In; i++)
    {
        x(1, i) = (i - 0.5) * dx; // cell-centered coordinates
        U.SetRow(i + nghost, InitialCondition(x(1, i)));
    }

    // set flux variable hold the memory
    matrix L(N, 8);
    matrix R(N, 8);
    matrix flux_L(1, 8);
    matrix flux_R(1, 8);
    matrix dflux(N, 8);
    matrix flux(N, 8);
    double dt;

    // open text
    ofstream myfile("anidata.txt");
    myfile << N_In << "\r\n";

    // main loop
    while(t < end_time)
//    for (int i = 0; i < largenumber; i++)
    {
        // set the boundary conditions
        BoundaryCondition(U);


        // estimate time-step from the CFL condition
        dt = ComputeTimestep(U);
        printf("t = %3.4f --> %3.4f\n", t, t + dt);

        // data reconstruction
        DataReconstruction_PLM(U, L, R);

        // update the face-centered variables by 0.5*dt
        for (int j = 2; j <= N - 1; j++)
        {
            flux_L = Conserved2FluxX(L.GetRow(j));
            flux_R = Conserved2FluxX(R.GetRow(j));
            dflux.SetRow(j, 0.5 * dt/dx * (flux_R-flux_L));
            L.SetRow(j, L.GetRow(j) - dflux.GetRow(j));
            R.SetRow(j, R.GetRow(j) - dflux.GetRow(j));
        }

        // compute fluxes
        // R[j-1] is the LEFT state at the j+1/2 inteface
#       pragma omp parallel for
        for (int j = nghost + 1; j <= N - nghost+1; j++)
        {
           // flux.SetRow(j, Roe(R.GetRow(j-1), L.GetRow(j), U(j, 6))); 
            flux.SetRow(j, HLLD(R.GetRow(j-1), L.GetRow(j), U(j, 6)));
        }
#       pragma omp parallel for
        // update the volume-averaged input variables by dt
        for (int j = nghost + 1; j <= N - nghost; j++)
        {
            U.SetRow(j, U.GetRow(j) - dt / dx * (flux.GetRow(j+1) - flux.GetRow(j)));
        }

        // update time
        t = t + dt;

        // write text
        if (myfile.is_open())
        {
            for (int j = nghost + 1; j <= N - nghost; j++)
            {
                myfile << t << "\t";
                myfile << x(1, j - 2) << "\t";
                myfile << U(j, 1) << "\t";
                myfile << U(j, 2)/U(j, 1) << "\t";
                myfile << ComputePressure(U.GetRow(j)) << "\t";
                myfile << "\r\n";
            }
        }
        else
            cout << "Unable to open file";
    
    }
    myfile.close();
    return 0;
}
