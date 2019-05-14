#include <iostream>
#include <stdlib.h>
#include <math.h>
// #include <cmath>
#include <stdio.h>

#include "ConstPara.h"
#include "matrix.h"
#include "physvar.h"

#include "Subfunc.cpp"
#include "InitBdaryCond.cpp"
#include "DataReconstruct.cpp"
#include "ComputeFlux.cpp"



main(int argc, char *argv[])
{
    physvar U;
    physvar L, R;
    matrix x(1, N_In);
    physvar flux_L;
    physvar flux_R;
    matrix dflux;
    physvar flux;
    x = InitialSpace(x);
    U = InitialCondition(x);

    double dt;

    for (int i = 0; i < largenumber; i++)
    {
        // set the boundary conditions
        U = BoundaryCondition(U);

        // estimate time-step from the CFL condition
        dt = 0.002;

        // data reconstruction
        DataReconstruction_PLM(U, L, R);

        // update the face-centered variables by 0.5*dt
        flux_L = Conserved2Flux(L);
        flux_R = Conserved2Flux(R);

        dflux = ComputedfluxSession(flux_L.d, flux_R.d, dt);
        L.d = L.d - dflux;
        R.d = R.d - dflux;
        dflux = ComputedfluxSession(flux_L.u, flux_R.u, dt);
        L.u = L.u - dflux;
        R.u = R.u - dflux;
        dflux = ComputedfluxSession(flux_L.E, flux_R.E, dt);
        L.E = L.E - dflux;
        R.E = R.E - dflux;
        
        // compute fluxes
        flux = Roe(L, R);

        // update the volume-averaged input variables by dt
        U = ComputeNext(flux, U, dt);
    }
    // U.d.display();
    // U.u.display();
    // U.E.display();
    return 0;
}