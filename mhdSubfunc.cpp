#ifndef __SUBFUNC_CPP__
#define __SUBFUNC_CPP__

#include "ConstPara.h"

// ##########################
// Compute pressure
// ##########################
double ComputePressure(matrix U)
{
    double P;
    P = (mygamma - 1.0)*( U(1,5) - 0.5*(U(1,2)*U(1,2) + U(1,3)*U(1,3) + U(1,4)*U(1,4))/U(1,1) - 0.5*(U(1,6)*U(1,6) + U(1,7)*U(1,7) + U(1,8)*U(1,8)) );
    if (P <= 0)
    {
        printf("negative pressure in compute pressure !!\t%f\n", P);
        exit(1);
    }
    return P;
}

double ComputePressureStar(const matrix &U)
{
    double PStar;
    PStar = (mygamma - 1.0)*( U(1,5) - 0.5*(U(1,2)*U(1,2) + U(1,3)*U(1,3) + U(1,4)*U(1,4))/U(1,1) );
    return PStar;
}

double ComputeEnergy(matrix W)
{
    double E;
    E = W(1,5)/(mygamma - 1.0) + ( 0.5*(W(1,2)*W(1,2) + W(1,3)*W(1,3) + W(1,4)*W(1,4))*W(1,1) + 0.5*(W(1,6)*W(1,6) + W(1,7)*W(1,7) + W(1,8)*W(1,8)) );
    return E;
}

// ##########################
// Convert conserved variables to primitive variables
// ##########################
matrix Conserved2Primitive(const matrix &U)
{
    matrix W(1, 8);
    W(1, 1) = U(1, 1);
    W(1, 2) = U(1, 2) / U(1, 1);
    W(1, 3) = U(1, 3) / U(1, 1);
    W(1, 4) = U(1, 4) / U(1, 1);
    W(1, 5) = ComputePressure(U);
    W(1, 6) = U(1, 6);
    W(1, 7) = U(1, 7);
    W(1, 8) = U(1, 8);
    return W;
}

// ##########################
// Convert primitive variables to conserved variables
// ##########################
matrix Primitive2Conserved(const matrix &W)
{
    matrix U(1, 8);
    U(1, 1) = W(1, 1);
    U(1, 2) = W(1, 1)*W(1, 2);
    U(1, 3) = W(1, 1)*W(1, 3);
    U(1, 4) = W(1, 1)*W(1, 4);
    U(1, 5) = ComputeEnergy(W);
    U(1, 6) = W(1, 6);
    U(1, 7) = W(1, 7);
    U(1, 8) = W(1, 8);
    return U;
}

// ##########################
// Define initial condition
// ##########################
matrix InitialCondition(double x)
{ // Sod shock tube
    matrix temp(1, 8);
    matrix U(1, 8);
    // define initial primiive variables
    if ( x < 0.5*L )
    {
        temp(1, 1) = 1.0;  // density
        temp(1, 2) = 0.0;  // velocity x
        temp(1, 3) = 0.0;  // velocity y
        temp(1, 4) = 0.0;  // velocity z
        temp(1, 5) = 1.0;  // pressure        
        temp(1, 6) = 0.75; // magnetic field x
        temp(1, 7) = 1.0;  // magnetic field y
        temp(1, 8) = 0.0;  // magnetic field z 
    }
    else
    {
        temp(1, 1) =  0.125; // density
        temp(1, 2) =  0.0;   // velocity x
        temp(1, 3) =  0.0;   // velocity y
        temp(1, 4) =  0.0;   // velocity z
        temp(1, 5) =  0.1;   // pressure        
        temp(1, 6) =  0.75;  // magnetic field x
        temp(1, 7) = -1.0;   // magnetic field y
        temp(1, 8) =  0.0;   // magnetic field z 
    }
    
    // conserved variables [1/2/3/4/5/6/7/8] <--> [density/momentum x/momentum y/momentum z/energy/Bx/By/Bz]
    U = Primitive2Conserved(temp); 
    return U;
}
/*
matrix DInitialCondition(double x, double y, double z )
{ // Sod shock tube
    matrix temp(1, 8);
    // define initial primiive variables
    if ( x < 0.5*L && y < 0.5*L && z < 0.5*L )
    {
        temp(1, 1) = 1.0;  // density
        temp(1, 2) = 0.0;  // velocity x
        temp(1, 3) = 0.0;  // velocity y
        temp(1, 4) = 0.0;  // velocity z
        temp(1, 5) = 1.0;  // pressure        
        temp(1, 6) = 0.75; // magnetic field x
        temp(1, 7) = 1.0;  // magnetic field y
        temp(1, 8) = 0.0;  // magnetic field z 
    }
    else
    {
        temp(1, 1) =  0.125; // density
        temp(1, 2) =  0.0;   // velocity x
        temp(1, 3) =  0.0;   // velocity y
        temp(1, 4) =  0.0;   // velocity z
        temp(1, 5) =  0.1;   // pressure        
        temp(1, 6) =  0.75;  // magnetic field x
        temp(1, 7) = -1.0;   // magnetic field y
        temp(1, 8) =  0.0;   // magnetic field z 
    }
    // conserved variables [1/2/3/4/5/6/7/8] <--> [density/momentum x/momentum y/momentum z/energy/Bx/By/Bz]
    temp = Primitive2Conserved(temp); 
    return temp;
}
*/

// ##########################
// Define boundary condition by setting ghost zones
// ##########################
void BoundaryCondition(matrix &U)
{ // outflow
    for (int i = 1; i <= nghost; i++)
    {
        for (int j = 1; j <= 8; j++)
        {
            U(i, j) = U(nghost + 1, j);
            U(N - i + 1, j) = U(N - nghost, j);
        }
    }
}

  
// ##########################
// Compute time-step by the CFL condition (directional unsplitting)
// ##########################
double ComputeTimestep(matrix &U)
{
    double dt_cfl = 0.0;
    double max_info_speed = 0.0;
    double P, a, u, v, w;
    double B2, cf;
    for (int i = 1; i <= N; i++)
    {
        P  = ComputePressure(U.GetRow(i));
        B2 = U(i, 6)*U(i, 6) + U(i, 7)*U(i, 7) + U(i, 8)*U(i, 8);
        a  = pow((mygamma * P / U(i, 1)), 0.5);
        cf = pow(0.5*(a*a + B2/U(i, 1)) + 0.5*pow( (a*a + B2/U(i, 1))*(a*a + B2/U(i, 1)) - 4.0*a*a*U(i, 6)*U(i, 6)/U(i ,1), 0.5) ,0.5);
        u  = abs(U(i, 2) / U(i, 1));
        v  = abs(U(i, 3) / U(i, 1));
        w  = abs(U(i, 4) / U(i, 1));
        if (max_info_speed < u + cf)
        {
            max_info_speed = u + cf;
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
    matrix slope_L(1, 8);
    matrix slope_R(1, 8);
    matrix slope_LR(1, 8);
    matrix slope_limited(1, 8);
    slope_L = C - L;
    slope_R = R - C;

    // apply the van-Leer limiter
    slope_LR = dotprod(slope_L, slope_R);
    for (int i = 1; i <= 8; i++)
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
// Piecewise-linear data reconstruction (PLM)
// ##########################
void DataReconstruction_PLM(matrix &U, matrix &L, matrix &R)
{
    // allocate memory
    matrix W(N, 8);
    matrix slope(N, 8);

    // conserved variables-- > primitive variables
    for (int j = 1; j <= N; j++)
        W.SetRow(j, Conserved2Primitive(U.GetRow(j)));

    // compute the left and right states of each cell
    for (int j = 2; j <= N - 1; j++)
        slope.SetRow(j, ComputeLimitedSlope(W.GetRow(j - 1), W.GetRow(j), W.GetRow(j + 1)));

    matrix slope_limited(1, 8);
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
    matrix W(N, 8);
    matrix slope(N, 8);
    matrix acc(N, 8);

    // conserved variables-- > primitive variables
    for (int j = 1; j <= N; j++)
        W.SetRow(j, Conserved2Primitive(U.GetRow(j)));

    // compute the left and right states of each cell
    for (int j = 2; j <= N - 1; j++)
        slope.SetRow(j, ComputeLimitedSlope(W.GetRow(j - 1), W.GetRow(j), W.GetRow(j + 1)));

    // compute the left and right states of each cell
    for (int j = 3; j <= N - 2; j++)
        acc.SetRow(j, ComputeLimitedSlope(slope.GetRow(j - 1), slope.GetRow(j), slope.GetRow(j + 1)));

    matrix slope_limited(1, 8);
    matrix acc_limited(1, 8);
    for (int j = 2; j <= N - 1; j++)
    {
        // get the face-centered variables
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
matrix Conserved2FluxX(matrix U)
{
    matrix fluxX(1, 8);
    double PStar = ComputePressureStar(U);
    double u = U(1, 2) / U(1, 1);
    double v = U(1, 3) / U(1, 1);
    double w = U(1, 4) / U(1, 1);
    fluxX(1, 1) = U(1, 2);
    fluxX(1, 2) = u*U(1, 2) - U(1, 6)*U(1, 6) + PStar;
    fluxX(1, 3) = u*U(1, 3) - U(1, 6)*U(1, 7);
    fluxX(1, 4) = u*U(1, 4) - U(1, 6)*U(1, 8);
    fluxX(1, 5) = u*( U(1, 5) + PStar ) - U(1, 6)*(u*U(1, 6) + v*U(1, 7) + w*U(1, 8));
    fluxX(1, 6) = 0;
    fluxX(1, 7) = u*U(1, 7) - v*U(1, 6);
    fluxX(1, 8) = u*U(1, 8) - w*U(1, 6);
    return fluxX;
}

matrix Conserved2FluxY(matrix U)
{
    matrix fluxY(1, 8);
    double PStar = ComputePressureStar(U);
    double u = U(1, 2) / U(1, 1);
    double v = U(1, 3) / U(1, 1);
    double w = U(1, 4) / U(1, 1);
    fluxY(1, 1) = U(1, 3);
    fluxY(1, 2) = v*U(1, 2) - U(1, 7)*U(1, 6);
    fluxY(1, 3) = v*U(1, 3) - U(1, 7)*U(1, 7) + PStar;
    fluxY(1, 4) = v*U(1, 4) - U(1, 7)*U(1, 8);
    fluxY(1, 5) = v*( U(1, 5) + PStar ) - U(1, 7)*(u*U(1, 6) + v*U(1, 7) + w*U(1, 8));
    fluxY(1, 6) = v*U(1, 6) - u*U(1, 7);
    fluxY(1, 7) = 0;
    fluxY(1, 8) = v*U(1, 8) - w*U(1, 7);
    return fluxY;
}

matrix Conserved2FluxZ(matrix U)
{
    matrix fluxZ(1, 8);
    double PStar = ComputePressureStar(U);
    double u = U(1, 2) / U(1, 1);
    double v = U(1, 3) / U(1, 1);
    double w = U(1, 4) / U(1, 1);
    fluxZ(1, 1) = U(1, 4);
    fluxZ(1, 2) = w*U(1, 2) - U(1, 8)*U(1, 6);
    fluxZ(1, 3) = w*U(1, 3) - U(1, 8)*U(1, 7);
    fluxZ(1, 4) = w*U(1, 4) - U(1, 8)*U(1, 8) + PStar;
    fluxZ(1, 5) = w*( U(1, 5) + PStar ) - U(1, 8)*(u*U(1, 6) + v*U(1, 7) + w*U(1, 8));
    fluxZ(1, 6) = w*U(1, 6) - u*U(1, 8);
    fluxZ(1, 7) = w*U(1, 7) - v*U(1, 8);
    fluxZ(1, 8) = 0;
    return fluxZ;
}

// ##########################
// Roe's Riemann solver for MHD
// ##########################
matrix RoeX(matrix L, matrix R, double Bx)
{
    // compute the enthalpy of the left and right states: H = (E+P)/rho
    matrix fluxX(1, 8);
    double P_L = ComputePressure(L);
    double P_R = ComputePressure(R);
    double PStar_L = ComputePressureStar(L);
    double PStar_R = ComputePressureStar(R);
    double H_L = (L(1, 5) + PStar_L) / L(1, 1);
    double H_R = (R(1, 5) + PStar_R) / R(1, 1);
    
    // compute Roe average values
    double rhoL_sqrt = pow(L(1, 1), 0.5);
    double rhoR_sqrt = pow(R(1, 1), 0.5);
    double rho = rhoL_sqrt*rhoR_sqrt;
    double u  = (L(1, 2) / rhoL_sqrt + R(1, 2) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double v  = (L(1, 3) / rhoL_sqrt + R(1, 3) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double w  = (L(1, 4) / rhoL_sqrt + R(1, 4) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);    
    double By = (L(1, 7) * rhoL_sqrt + R(1, 7) * rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double Bz = (L(1, 8) * rhoL_sqrt + R(1, 8) * rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double H  = (rhoL_sqrt * H_L + rhoR_sqrt * H_R) / (rhoL_sqrt + rhoR_sqrt);
    double V2 = u*u + v*v + w*w;
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double X  = 0.5*(pow(R(1, 7) - L(1, 7), 2.0) + pow(R(1, 8) - L(1, 8), 2.0))/pow(rhoL_sqrt + rhoR_sqrt, 2.0); 
    double P  = (rho/mygamma)*((mygamma-1.0)*(H - 0.5*V2 - B2/rho ) - (mygamma-2.0)*X);

    // check negative pressure
    if ( P < 0.0 )
    {
        printf("negative pressure!\t%f\n", P);
        exit(1);
    }
    // calculate speed 
    double a  = pow( mygamma*P/rho ,0.5);  //sound speed
    double ca = pow( Bx*Bx/rho, 0.5);      //Alfven speed
    double cf = pow(0.5*(a*a + B2/rho) + 0.5*pow( (a*a + B2/rho)*(a*a + B2/rho) - 4.0*a*a*Bx*Bx/rho, 0.5) ,0.5); //fast magnetosonic speed
    double cs = pow(0.5*(a*a + B2/rho) - 0.5*pow( (a*a + B2/rho)*(a*a + B2/rho) - 4.0*a*a*Bx*Bx/rho, 0.5) ,0.5); //slow magnetosonic speed
    
    // compute some parameters
    double alphaf = pow((a*a - cs*cs)/(cf*cf - cs*cs), 0.5);
    double alphas = pow((cf*cf - a*a)/(cf*cf - cs*cs), 0.5); 
    double betay  = By/abs( pow(By*By + Bz+Bz, 0.5) );
    double betaz  = Bz/abs( pow(By*By + Bz+Bz, 0.5) );
    double rho_sqrt = pow(rho, 0.5);
    double S;    //S = sign(Bx)
    if(Bx > 0.0)
    {
        S = 1.0;
    }
    else if(Bx < 0.0)
    {
        S = -1.0;
    }
    else
    {
        S = 0.0;
    }

    // compute the amplitudes of different characteristic waves
    matrix dU(1, 8);
    matrix amp(1, 7);

    dU = R - L;
    
    amp(1, 1) = 0.5*(alphaf*(X*dU(1, 5) + (P_R - P_L)) + rho*alphas*cs*S*(betay*dU(1, 3) + betaz*dU(1, 4)) - rho*alphaf*cf*dU(1, 2) + rho_sqrt*alphas*a*(betay*dU(1, 7) + betaz*dU(1,8))); 
    amp(1, 2) = 0.5*(betay*dU(1, 4) - betaz*dU(1 ,3) + S*(betay*dU(1, 8)- betaz*dU(1, 7))/rho_sqrt); 
    amp(1, 3) = 0.5*(alphas*(X*dU(1, 5) + (P_R - P_L)) - rho*alphaf*cf*S*(betay*dU(1, 3) + betaz*dU(1, 4)) - rho*alphas*cs*dU(1, 2) + rho_sqrt*alphaf*a*(betay*dU(1, 7) + betaz*dU(1,8)));
    amp(1, 4) = (a*a - X)*dU(1, 1) - (P_R - P_L);
    amp(1, 5) = 0.5*(alphas*(X*dU(1, 5) + (P_R - P_L)) + rho*alphaf*cf*S*(betay*dU(1, 3) + betaz*dU(1, 4)) + rho*alphas*cs*dU(1, 2) + rho_sqrt*alphaf*a*(betay*dU(1, 7) + betaz*dU(1,8)));
    amp(1, 6) = 0.5*(betaz*dU(1, 3) - betay*dU(1 ,4) + S*(betay*dU(1, 8)- betaz*dU(1, 7))/rho_sqrt);
    amp(1, 7) = 0.5*(alphaf*(X*dU(1, 5) + (P_R - P_L)) - rho*alphas*cs*S*(betay*dU(1, 3) + betaz*dU(1, 4)) + rho*alphaf*cf*dU(1, 2) + rho_sqrt*alphas*a*(betay*dU(1, 7) + betaz*dU(1,8)));
    // compute the eigenvalues and right eigenvector matrix
    matrix EigVal(1, 7);
    matrix EigVec(7, 8);

    EigVal(1, 1) = abs(u - cf);
    EigVal(1, 2) = abs(u - ca);
    EigVal(1, 3) = abs(u - cs);
    EigVal(1, 4) = abs(u);
    EigVal(1, 5) = abs(u + cs);
    EigVal(1, 6) = abs(u + ca);
    EigVal(1, 7) = abs(u + cf);
    
    //R_u-cf 
    EigVec(1, 1) = alphaf/(a*a);
    EigVec(1, 2) = alphaf*(u - cf)/(a*a);
    EigVec(1, 3) = (alphaf*v + alphas*cs*betay*S)/(a*a);
    EigVec(1, 4) = (alphaf*w + alphas*cs*betaz*S)/(a*a);
    EigVec(1, 5) = ( alphaf*(H - B2/rho - u*cf) + alphas*cs*S*(v*betay + w*betaz) - alphas*a*abs(pow(By*By + Bz+Bz, 0.5))/rho_sqrt )/(a*a);
    EigVec(1, 6) = 0.0;
    EigVec(1, 7) = alphas*a*betay/(rho_sqrt*a*a);
    EigVec(1, 8) = alphas*a*betaz/(rho_sqrt*a*a);
    //R_u-ca
    EigVec(2, 1) = 0.0;
    EigVec(2, 2) = 0.0;
    EigVec(2, 3) = -rho*betaz;
    EigVec(2, 4) =  rho*betay;
    EigVec(2, 5) = -rho*(v*betaz - w*betay);
    EigVec(2, 6) = 0.0;
    EigVec(2, 7) = -S*rho_sqrt*betaz;
    EigVec(2, 8) =  S*rho_sqrt*betay;
    //R_u-cs
    EigVec(3, 1) = alphas/(a*a);
    EigVec(3, 2) = alphas*(u - cs)/(a*a);
    EigVec(3, 3) = (alphas*v - alphaf*cf*betay*S)/(a*a);
    EigVec(3, 4) = (alphas*w - alphaf*cf*betaz*S)/(a*a);
    EigVec(3, 5) = ( alphas*(H - B2/rho - u*cs) - alphaf*cf*S*(v*betay + w*betaz) - alphaf*a*abs(pow(By*By + Bz+Bz, 0.5))/rho_sqrt )/(a*a);
    EigVec(3, 6) = 0.0;
    EigVec(3, 7) = -alphaf*a*betay/(rho_sqrt*a*a);
    EigVec(3, 8) = -alphaf*a*betaz/(rho_sqrt*a*a);
    //R_u
    EigVec(4, 1) = 1.0/(a*a); 
    EigVec(4, 2) = u/(a*a);
    EigVec(4, 3) = v/(a*a);
    EigVec(4, 4) = w/(a*a);
    EigVec(4, 5) = (0.5*V2 + X*(mygamma-2.0)/(mygamma-1.0))/(a*a); 
    EigVec(4, 6) = 0.0;
    EigVec(4, 7) = 0.0;
    EigVec(4, 8) = 0.0;
    //R_u+cs
    EigVec(5, 1) = alphas/(a*a); 
    EigVec(5, 2) = alphas*(u + cs)/(a*a);
    EigVec(5, 3) = (alphas*v + alphaf*cf*betay*S)/(a*a);
    EigVec(5, 4) = (alphas*w + alphaf*cf*betaz*S)/(a*a);
    EigVec(5, 5) = ( alphas*(H - B2/rho + u*cs) + alphaf*cf*S*(v*betay + w*betaz) - alphaf*a*abs(pow(By*By + Bz+Bz, 0.5))/rho_sqrt )/(a*a);
    EigVec(5, 6) = 0.0;
    EigVec(5, 7) = -alphaf*a*betay/(rho_sqrt*a*a);
    EigVec(5, 8) = -alphaf*a*betaz/(rho_sqrt*a*a);
    //R_u+ca
    EigVec(6, 1) = 0.0;
    EigVec(6, 2) = 0.0;
    EigVec(6, 3) =  rho*betaz;
    EigVec(6, 4) = -rho*betay;
    EigVec(6, 5) =  rho*(v*betaz - w*betay);
    EigVec(6, 6) = 0.0;
    EigVec(6, 7) = -S*rho_sqrt*betaz;
    EigVec(6, 8) =  S*rho_sqrt*betay;
    //R_u+cf
    EigVec(7, 1) = alphaf/(a*a);
    EigVec(7, 2) = alphaf*(u + cf)/(a*a);
    EigVec(7, 3) = (alphaf*v - alphas*cs*betay*S)/(a*a);
    EigVec(7, 4) = (alphaf*w - alphas*cs*betaz*S)/(a*a);
    EigVec(7, 5) = ( alphaf*(H - B2/rho + u*cf) - alphas*cs*S*(v*betay + w*betaz) - alphas*a*abs(pow(By*By + Bz+Bz, 0.5))/rho_sqrt )/(a*a);
    EigVec(7, 6) = 0.0;
    EigVec(7, 7) = alphas*a*betay/(rho_sqrt*a*a);
    EigVec(7, 8) = alphas*a*betaz/(rho_sqrt*a*a);
    
    // compute the fluxes of the left and right states
    matrix flux_L(1, 8);
    matrix flux_R(1, 8);

    flux_L = Conserved2FluxX(L);
    flux_R = Conserved2FluxX(R);

    // compute the Roe flux
    for(int i = 1; i <= 8; i++)
    {
        for(int j = 1; j <= 8; j++)
        {
            fluxX(1, i) = 0.5*(flux_L(1, i) + flux_R(1, i)) - 0.5*amp(1, j)*EigVal(1, j)*EigVec(j, i);
        }
    } 
    return fluxX;
}

matrix RoeY(matrix L, matrix R, double By)
{
    // compute the enthalpy of the left and right states: H = (E+P*)/rho
    matrix fluxY(1, 8);
    double P_L = ComputePressure(L);
    double P_R = ComputePressure(R);
    double PStar_L = ComputePressureStar(L);
    double PStar_R = ComputePressureStar(R);
    double H_L = (L(1, 5) + PStar_L) / L(1, 1);
    double H_R = (R(1, 5) + PStar_R) / R(1, 1);
    
    // compute Roe average values
    double rhoL_sqrt = pow(L(1, 1), 0.5);
    double rhoR_sqrt = pow(R(1, 1), 0.5);
    double rho = rhoL_sqrt*rhoR_sqrt;
    double u  = (L(1, 2) / rhoL_sqrt + R(1, 2) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double v  = (L(1, 3) / rhoL_sqrt + R(1, 3) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double w  = (L(1, 4) / rhoL_sqrt + R(1, 4) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);    
    double Bx = (L(1, 6) * rhoL_sqrt + R(1, 6) * rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double Bz = (L(1, 8) * rhoL_sqrt + R(1, 8) * rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double H  = (rhoL_sqrt * H_L + rhoR_sqrt * H_R) / (rhoL_sqrt + rhoR_sqrt);
    double V2 = u*u + v*v + w*w;
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double X  = 0.5*(pow(R(1, 6) - L(1, 6), 2.0) + pow(R(1, 8) - L(1, 8), 2.0))/pow(rhoL_sqrt + rhoR_sqrt, 2.0); 
    double P  = (rho/mygamma)*((mygamma-1.0)*(H - 0.5*V2 - B2/rho ) - (mygamma-2.0)*X);

    // check negative pressure
    if ( P < 0.0 )
    {
        printf("negative pressure!\t%f\n", H - 0.5 * V2);
        exit(1);
    }
    // calculate speed 
    double a  = pow( mygamma*P/rho ,0.5);  //sound speed
    double ca = pow( By*By/rho, 0.5);      //Alfven speed
    double cf = pow(0.5*(a*a + B2/rho) + 0.5*pow( (a*a + B2/rho)*(a*a + B2/rho) - 4.0*a*a*By*By/rho, 0.5) ,0.5); //fast magnetosonic speed
    double cs = pow(0.5*(a*a + B2/rho) - 0.5*pow( (a*a + B2/rho)*(a*a + B2/rho) - 4.0*a*a*By*By/rho, 0.5) ,0.5); //slow magnetosonic speed
    
    // compute some parameters
    double alphaf = pow((a*a - cs*cs)/(cf*cf - cs*cs), 0.5);
    double alphas = pow((cf*cf - a*a)/(cf*cf - cs*cs), 0.5); 
    double betax  = Bx/abs( pow(Bx*Bx + Bz+Bz, 0.5) );
    double betaz  = Bz/abs( pow(Bx*Bx + Bz+Bz, 0.5) );
    double rho_sqrt = pow(rho, 0.5);
    double S;    //S = sign(By)
    if(By > 0.0)
    {
        S = 1.0;
    }
    else if(By < 0.0)
    {
        S = -1.0;
    }
    else
    {
        S = 0.0;
    }

    // compute the amplitudes of different characteristic waves
    matrix dU(1, 8);
    matrix amp(1,7);

    dU = R - L;
    
    amp(1, 1) = 0.5*(alphaf*(X*dU(1, 5) + (P_R - P_L)) + rho*alphas*cs*S*(betaz*dU(1, 4) + betax*dU(1, 2)) - rho*alphaf*cf*dU(1, 3) + rho_sqrt*alphas*a*(betaz*dU(1, 8) + betax*dU(1, 6))); 
    amp(1, 2) = 0.5*(betaz*dU(1, 2) - betax*dU(1 ,4) + S*(betaz*dU(1, 6)- betax*dU(1, 8))/rho_sqrt); 
    amp(1, 3) = 0.5*(alphas*(X*dU(1, 5) + (P_R - P_L)) - rho*alphaf*cf*S*(betaz*dU(1, 4) + betax*dU(1, 2)) - rho*alphas*cs*dU(1, 3) + rho_sqrt*alphaf*a*(betaz*dU(1, 8) + betax*dU(1, 6)));
    amp(1, 4) = (a*a - X)*dU(1, 1) - (P_R - P_L);
    amp(1, 5) = 0.5*(alphas*(X*dU(1, 5) + (P_R - P_L)) + rho*alphaf*cf*S*(betaz*dU(1, 4) + betax*dU(1, 2)) + rho*alphas*cs*dU(1, 3) + rho_sqrt*alphaf*a*(betaz*dU(1, 8) + betax*dU(1, 6)));
    amp(1, 6) = 0.5*(betax*dU(1, 4) - betaz*dU(1 ,2) + S*(betaz*dU(1, 6)- betax*dU(1, 8))/rho_sqrt);
    amp(1, 7) = 0.5*(alphaf*(X*dU(1, 5) + (P_R - P_L)) - rho*alphas*cs*S*(betaz*dU(1, 4) + betax*dU(1, 2)) + rho*alphaf*cf*dU(1, 3) + rho_sqrt*alphas*a*(betaz*dU(1, 8) + betax*dU(1, 6)));
    // compute the eigenvalues and right eigenvector matrix R
    matrix EigVal(1, 7);
    matrix EigVec(7, 8);

    EigVal(1, 1) = abs(v - cf);
    EigVal(1, 2) = abs(v - ca);
    EigVal(1, 3) = abs(v - cs);
    EigVal(1, 4) = abs(v);
    EigVal(1, 5) = abs(v + cs);
    EigVal(1, 6) = abs(v + ca);
    EigVal(1, 7) = abs(v + cf);
    
    //R_v-cf 
    EigVec(1, 1) = alphaf/(a*a);
    EigVec(1, 2) = (alphaf*u + alphas*cs*betax*S)/(a*a);
    EigVec(1, 3) = alphaf*(v - cf)/(a*a);
    EigVec(1, 4) = (alphaf*w + alphas*cs*betaz*S)/(a*a);
    EigVec(1, 5) = ( alphaf*(H - B2/rho - v*cf) + alphas*cs*S*(u*betax + w*betaz) - alphas*a*abs(pow(Bx*Bx + Bz+Bz, 0.5))/rho_sqrt )/(a*a);
    EigVec(1, 6) = alphas*a*betax/(rho_sqrt*a*a);
    EigVec(1, 7) = 0.0;
    EigVec(1, 8) = alphas*a*betaz/(rho_sqrt*a*a);
    //R_v-ca
    EigVec(2, 1) = 0.0;
    EigVec(2, 2) =  rho*betaz;
    EigVec(2, 3) = 0.0;
    EigVec(2, 4) = -rho*betax;
    EigVec(2, 5) = -rho*(w*betax - u*betaz);
    EigVec(2, 6) =  S*rho_sqrt*betaz;
    EigVec(2, 7) = 0.0;
    EigVec(2, 8) = -S*rho_sqrt*betax;
    //R_v-cs
    EigVec(3, 1) = alphas/(a*a);
    EigVec(3, 2) = (alphas*u - alphaf*cf*betax*S)/(a*a);
    EigVec(3, 3) = alphas*(v - cs)/(a*a);
    EigVec(3, 4) = (alphas*w - alphaf*cf*betaz*S)/(a*a);
    EigVec(3, 5) = ( alphas*(H - B2/rho - v*cs) - alphaf*cf*S*(w*betaz + u*betax) - alphaf*a*abs(pow(Bx*Bx + Bz+Bz, 0.5))/rho_sqrt )/(a*a);
    EigVec(3, 6) = -alphaf*a*betax/(rho_sqrt*a*a);
    EigVec(3, 7) = 0.0;
    EigVec(3, 8) = -alphaf*a*betaz/(rho_sqrt*a*a);
    //R_v
    EigVec(4, 1) = 1.0/(a*a); 
    EigVec(4, 2) = u/(a*a);
    EigVec(4, 3) = v/(a*a);
    EigVec(4, 4) = w/(a*a);
    EigVec(4, 5) = (0.5*V2 + X*(mygamma-2.0)/(mygamma-1.0))/(a*a); 
    EigVec(4, 6) = 0.0;
    EigVec(4, 7) = 0.0;
    EigVec(4, 8) = 0.0;
    //R_v+cs
    EigVec(5, 1) = alphas/(a*a);
    EigVec(5, 2) = (alphas*u + alphaf*cf*betax*S)/(a*a);
    EigVec(5, 3) = alphas*(v + cs)/(a*a);
    EigVec(5, 4) = (alphas*w + alphaf*cf*betaz*S)/(a*a);
    EigVec(5, 5) = ( alphas*(H - B2/rho + v*cs) + alphaf*cf*S*(w*betaz + u*betax) - alphaf*a*abs(pow(Bx*Bx + Bz+Bz, 0.5))/rho_sqrt )/(a*a);
    EigVec(5, 6) = -alphaf*a*betax/(rho_sqrt*a*a);
    EigVec(5, 7) = 0.0;
    EigVec(5, 8) = -alphaf*a*betaz/(rho_sqrt*a*a);
    //R_v+ca
    EigVec(6, 1) = 0.0;
    EigVec(6, 2) = -rho*betaz;
    EigVec(6, 3) = 0.0;
    EigVec(6, 4) =  rho*betax;
    EigVec(6, 5) =  rho*(w*betax - u*betaz);
    EigVec(6, 6) =  S*rho_sqrt*betaz;
    EigVec(6, 7) = 0.0;
    EigVec(6, 8) = -S*rho_sqrt*betax;
    //R_v+cf
    EigVec(7, 1) = alphaf/(a*a);
    EigVec(7, 2) = (alphaf*u - alphas*cs*betax*S)/(a*a);
    EigVec(7, 3) = alphaf*(v + cf)/(a*a);
    EigVec(7, 4) = (alphaf*w - alphas*cs*betaz*S)/(a*a);
    EigVec(7, 5) = ( alphaf*(H - B2/rho + v*cf) - alphas*cs*S*(u*betax + w*betaz) - alphas*a*abs(pow(Bx*Bx + Bz+Bz, 0.5))/rho_sqrt )/(a*a);
    EigVec(7, 6) = alphas*a*betax/(rho_sqrt*a*a);
    EigVec(7, 7) = 0.0;
    EigVec(7, 8) = alphas*a*betaz/(rho_sqrt*a*a);
    
    // compute the fluxes of the left and right states
    matrix flux_L(1, 8);
    matrix flux_R(1, 8);

    flux_L = Conserved2FluxY(L);
    flux_R = Conserved2FluxY(R);

    // compute the Roe flux
    for(int i = 1; i <= 8; i++)                                                                           
    {                                                                                                     
        for(int j = 1; j <= 8; j++)                                                                       
        { 
            fluxY(1, i) = 0.5*(flux_L(1, i) + flux_R(1, i)) - 0.5*amp(1, j)*EigVal(1, j)*EigVec(j, i);    
        }
    } 
    return fluxY;
}

matrix RoeZ(matrix L, matrix R, double Bz)
{
    // compute the enthalpy of the left and right states: H = (E+P*)/rho
    matrix fluxZ(1, 8);
    double P_L = ComputePressure(L);
    double P_R = ComputePressure(R);
    double PStar_L = ComputePressureStar(L);
    double PStar_R = ComputePressureStar(R);
    double H_L = (L(1, 5) + PStar_L) / L(1, 1);
    double H_R = (R(1, 5) + PStar_R) / R(1, 1);
    
    // compute Roe average values
    double rhoL_sqrt = pow(L(1, 1), 0.5);
    double rhoR_sqrt = pow(R(1, 1), 0.5);
    double rho = rhoL_sqrt*rhoR_sqrt;
    double u  = (L(1, 2) / rhoL_sqrt + R(1, 2) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double v  = (L(1, 3) / rhoL_sqrt + R(1, 3) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double w  = (L(1, 4) / rhoL_sqrt + R(1, 4) / rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);    
    double Bx = (L(1, 6) * rhoL_sqrt + R(1, 6) * rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double By = (L(1, 7) * rhoL_sqrt + R(1, 7) * rhoR_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double H  = (rhoL_sqrt * H_L + rhoR_sqrt * H_R) / (rhoL_sqrt + rhoR_sqrt);
    double V2 = u*u + v*v + w*w;
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double X  = 0.5*(pow(R(1, 6) - L(1, 6), 2.0) + pow(R(1, 7) - L(1, 7), 2.0))/pow(rhoL_sqrt + rhoR_sqrt, 2.0); 
    double P  = (rho/mygamma)*((mygamma-1.0)*(H - 0.5*V2 - B2/rho ) - (mygamma-2.0)*X);

    // check negative pressure
    if ( P < 0.0 )
    {
        printf("negative pressure!\t%f\n", H - 0.5 * V2);
        exit(1);
    }
    // calculate speed 
    double a  = pow( mygamma*P/rho ,0.5);  //sound speed
    double ca = pow( Bz*Bz/rho, 0.5);      //Alfven speed
    double cf = pow(0.5*(a*a + B2/rho) + 0.5*pow( (a*a + B2/rho)*(a*a + B2/rho) - 4.0*a*a*Bz*Bz/rho, 0.5) ,0.5); //fast magnetosonic speed
    double cs = pow(0.5*(a*a + B2/rho) - 0.5*pow( (a*a + B2/rho)*(a*a + B2/rho) - 4.0*a*a*Bz*Bz/rho, 0.5) ,0.5); //slow magnetosonic speed
    
    // compute some parameters
    double alphaf = pow((a*a - cs*cs)/(cf*cf - cs*cs), 0.5);
    double alphas = pow((cf*cf - a*a)/(cf*cf - cs*cs), 0.5); 
    double betax  = Bx/abs( pow(Bx*Bx + By+By, 0.5) );
    double betay  = By/abs( pow(Bx*Bx + By+By, 0.5) );
    double rho_sqrt = pow(rho, 0.5);
    double S;    //S = sign(Bz)
    if(Bz > 0.0)
    {
        S = 1.0;
    }
    else if(Bz < 0.0)
    {
        S = -1.0;
    }
    else
    {
        S = 0.0;
    }

    // compute the amplitudes of different characteristic waves
    matrix dU(1, 8);
    matrix amp(1,7);

    dU = R - L;
    
    amp(1, 1) = 0.5*(alphaf*(X*dU(1, 5) + (P_R - P_L)) + rho*alphas*cs*S*(betax*dU(1, 2) + betay*dU(1, 3)) - rho*alphaf*cf*dU(1, 4) + rho_sqrt*alphas*a*(betax*dU(1, 6) + betay*dU(1, 7))); 
    amp(1, 2) = 0.5*(betax*dU(1, 3) - betay*dU(1 ,2) + S*(betax*dU(1, 7)- betay*dU(1, 6))/rho_sqrt); 
    amp(1, 3) = 0.5*(alphas*(X*dU(1, 5) + (P_R - P_L)) - rho*alphaf*cf*S*(betax*dU(1, 2) + betay*dU(1, 3)) - rho*alphas*cs*dU(1, 4) + rho_sqrt*alphaf*a*(betax*dU(1, 6) + betay*dU(1, 7)));
    amp(1, 4) = (a*a - X)*dU(1, 1) - (P_R - P_L);
    amp(1, 5) = 0.5*(alphas*(X*dU(1, 5) + (P_R - P_L)) + rho*alphaf*cf*S*(betax*dU(1, 2) + betay*dU(1, 3)) + rho*alphas*cs*dU(1, 4) + rho_sqrt*alphaf*a*(betax*dU(1, 6) + betay*dU(1, 7)));
    amp(1, 6) = 0.5*(betay*dU(1, 2) - betax*dU(1 ,3) + S*(betax*dU(1, 7)- betay*dU(1, 6))/rho_sqrt);
    amp(1, 7) = 0.5*(alphaf*(X*dU(1, 5) + (P_R - P_L)) - rho*alphas*cs*S*(betax*dU(1, 2) + betay*dU(1, 3)) + rho*alphaf*cf*dU(1, 4) + rho_sqrt*alphas*a*(betax*dU(1, 6) + betay*dU(1, 7)));
    // compute the eigenvalues and right eigenvector matrix R
    matrix EigVal(1, 7);
    matrix EigVec(7, 8);

    EigVal(1, 1) = abs(w - cf);
    EigVal(1, 2) = abs(w - ca);
    EigVal(1, 3) = abs(w - cs);
    EigVal(1, 4) = abs(w);
    EigVal(1, 5) = abs(w + cs);
    EigVal(1, 6) = abs(w + ca);
    EigVal(1, 7) = abs(w + cf);
    
    //R_w-cf 
    EigVec(1, 1) = alphaf/(a*a);
    EigVec(1, 2) = (alphaf*u + alphas*cs*betax*S)/(a*a);
    EigVec(1, 3) = (alphaf*v + alphas*cs*betay*S)/(a*a);
    EigVec(1, 4) = alphaf*(w - cf)/(a*a);
    EigVec(1, 5) = ( alphaf*(H - B2/rho - w*cf) + alphas*cs*S*(v*betay + u*betax) - alphas*a*abs(pow(Bx*Bx + By+By, 0.5))/rho_sqrt )/(a*a);
    EigVec(1, 6) = alphas*a*betax/(rho_sqrt*a*a);
    EigVec(1, 7) = alphas*a*betay/(rho_sqrt*a*a);
    EigVec(1, 8) = 0.0;
    //R_w-ca
    EigVec(2, 1) = 0.0;
    EigVec(2, 2) = -rho*betay;
    EigVec(2, 3) =  rho*betax;
    EigVec(2, 4) = 0.0;
    EigVec(2, 5) = -rho*(u*betay - v*betax);
    EigVec(2, 6) = -S*rho_sqrt*betay;
    EigVec(2, 7) =  S*rho_sqrt*betax;
    EigVec(2, 8) = 0.0;
    //R_w-cs
    EigVec(3, 1) = alphas/(a*a);
    EigVec(3, 2) = (alphas*u - alphaf*cf*betax*S)/(a*a);
    EigVec(3, 3) = (alphas*v - alphaf*cf*betay*S)/(a*a);
    EigVec(3, 4) = alphas*(w - cs)/(a*a);
    EigVec(3, 5) = ( alphas*(H - B2/rho - w*cs) - alphaf*cf*S*(u*betax + v*betay) - alphaf*a*abs(pow(Bx*Bx + By+By, 0.5))/rho_sqrt )/(a*a);
    EigVec(3, 6) = -alphaf*a*betax/(rho_sqrt*a*a);
    EigVec(3, 7) = -alphaf*a*betay/(rho_sqrt*a*a);
    EigVec(3, 8) = 0.0;
    //R_w
    EigVec(4, 1) = 1.0/(a*a); 
    EigVec(4, 2) = u/(a*a);
    EigVec(4, 3) = v/(a*a);
    EigVec(4, 4) = w/(a*a);
    EigVec(4, 5) = (0.5*V2 + X*(mygamma-2.0)/(mygamma-1.0))/(a*a); 
    EigVec(4, 6) = 0.0;
    EigVec(4, 7) = 0.0;
    EigVec(4, 8) = 0.0;
    //R_w+cs
    EigVec(5, 1) = alphas/(a*a);
    EigVec(5, 2) = (alphas*u + alphaf*cf*betax*S)/(a*a);
    EigVec(5, 3) = (alphas*v + alphaf*cf*betay*S)/(a*a);
    EigVec(5, 4) = alphas*(w + cs)/(a*a);
    EigVec(5, 5) = ( alphas*(H - B2/rho + w*cs) + alphaf*cf*S*(u*betax + v*betay) - alphaf*a*abs(pow(Bx*Bx + By+By, 0.5))/rho_sqrt )/(a*a);
    EigVec(5, 6) = -alphaf*a*betax/(rho_sqrt*a*a);
    EigVec(5, 7) = -alphaf*a*betay/(rho_sqrt*a*a);
    EigVec(5, 8) = 0.0;
    //R_w+ca
    EigVec(6, 1) = 0.0;
    EigVec(6, 2) =  rho*betay;
    EigVec(6, 3) = -rho*betax;
    EigVec(6, 4) = 0.0;
    EigVec(6, 5) =  rho*(u*betay - v*betax);
    EigVec(6, 6) = -S*rho_sqrt*betay;
    EigVec(6, 7) =  S*rho_sqrt*betax;
    EigVec(6, 8) = 0.0;
    //R_w+cf
    EigVec(7, 1) = alphaf/(a*a);
    EigVec(7, 2) = (alphaf*u - alphas*cs*betax*S)/(a*a);
    EigVec(7, 3) = (alphaf*v - alphas*cs*betay*S)/(a*a);
    EigVec(7, 4) = alphaf*(w + cf)/(a*a);
    EigVec(7, 5) = ( alphaf*(H - B2/rho + w*cf) - alphas*cs*S*(v*betay + u*betax) - alphas*a*abs(pow(Bx*Bx + By+By, 0.5))/rho_sqrt )/(a*a);
    EigVec(7, 6) = alphas*a*betax/(rho_sqrt*a*a);
    EigVec(7, 7) = alphas*a*betay/(rho_sqrt*a*a);
    EigVec(7, 8) = 0.0;
    // compute the fluxes of the left and right states
    matrix flux_L(1, 8);
    matrix flux_R(1, 8);

    flux_L = Conserved2FluxZ(L);
    flux_R = Conserved2FluxZ(R);

    // compute the Roe flux
    for(int i = 1; i <= 8; i++)
    {
        for(int j = 1; j <= 8; j++)
        {
            fluxZ(1, i) = 0.5*(flux_L(1, i) + flux_R(1, i)) - 0.5*amp(1, j)*EigVal(1, j)*EigVec(j, i);
        }
    }
    return fluxZ;
}
#endif //__SUBFUNC_CPP__
