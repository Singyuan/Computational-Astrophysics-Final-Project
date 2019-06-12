#ifndef __SUBFUNC_CPP__
#define __SUBFUNC_CPP__

#include "ConstPara.h"
using namespace std;
// ##########################
// Compute pressure
// ##########################
double ComputePressure(const matrix &U)
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
    PStar = ComputePressure(U) + 0.5*(U(1,6)*U(1,6) + U(1,7)*U(1,7) + U(1,8)*U(1,8)); 
    return PStar;
}

double ComputeEnergy(const matrix &W)
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
        temp(1, 5) =  0.1;  // pressure        
        temp(1, 6) =  0.75;  // magnetic field x
        temp(1, 7) = -1.0;   // magnetic field y
        temp(1, 8) =  0.0;   // magnetic field z 
    }
    
    // conserved variables [1/2/3/4/5/6/7/8] <--> [density/momentum x/momentum y/momentum z/energy/Bx/By/Bz]
    U = Primitive2Conserved(temp); 
    return U;
}

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
matrix ComputeLimitedSlope(const matrix &L, const matrix &C, const matrix &R)
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
 
    // conserved variables-- > primitive variables
    for (int j = 1; j <= N; j++)
        W.SetRow(j, Conserved2Primitive(U.GetRow(j)));

    // compute the left and right states of each cell
    for (int j = 2; j <= N - 1; j++)
    {   
        // compute the left and right states of each cell
        matrix slope_limited(1, 8);
        slope_limited = ComputeLimitedSlope(W.GetRow(j - 1), W.GetRow(j), W.GetRow(j + 1));
    
        //  get the face-centered variables
        L.SetRow(j, W.GetRow(j) - 0.5 * slope_limited);
        R.SetRow(j, W.GetRow(j) + 0.5 * slope_limited);

        // ensure face-centered variables lie between nearby volume-averaged (~cell-centered) values
        for(int i = 0; i < 8; i++ )
        {
            L(j, i) = max( L(j, i), min( W(j-1, i), W(j, i) ) );
            L(j, i) = min( L(j, i), max( W(j-1, i), W(j, i) ) );
            R(j, i) = 2.0*W(j, i) - L(j, i);

            R(j, i) = max( R(j, i), min( W(j+1, i), W(j, i) ) );
            R(j, i) = min( R(j, i), max( W(j+1, i), W(j, i) ) );
            L(j, i) = 2.0*W(j, i) - R(j, i);
        } 
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
matrix Conserved2FluxX(const matrix &U)
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
    fluxX(1, 6) = 0.0;
    fluxX(1, 7) = u*U(1, 7) - v*U(1, 6);
    fluxX(1, 8) = u*U(1, 8) - w*U(1, 6);
    return fluxX;
}

matrix Conserved2FluxY(const matrix &U)
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
    fluxY(1, 7) = 0.0;
    fluxY(1, 8) = v*U(1, 8) - w*U(1, 7);
    return fluxY;
}

matrix Conserved2FluxZ(const matrix &U)
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
    fluxZ(1, 8) = 0.0;
    return fluxZ;
}

// ##########################
// Roe's Riemann solver for MHD
// ##########################
matrix Roe(const matrix &L, const matrix &R, double Bx)
{
    // compute the enthalpy of the left and right states: H = (E+P*)/rho
    matrix flux(1, 8);
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
    double By = (L(1, 7) * rhoR_sqrt + R(1, 7) * rhoL_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double Bz = (L(1, 8) * rhoR_sqrt + R(1, 8) * rhoL_sqrt) / (rhoL_sqrt + rhoR_sqrt);
    double H  = (rhoL_sqrt * H_L + rhoR_sqrt * H_R) / (rhoL_sqrt + rhoR_sqrt);
    double V2 = u*u + v*v + w*w;
    double B2 = Bx*Bx + By*By + Bz*Bz;
    double X  = 0.5*(pow(L(1, 7) - R(1, 7), 2.0) + pow(L(1, 8) - R(1, 8), 2.0))/pow(rhoL_sqrt + rhoR_sqrt, 2.0); 
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
    double betay  = By/pow(By*By + Bz*Bz, 0.5);
    double betaz  = Bz/pow(By*By + Bz*Bz, 0.5);
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
    matrix dW(1, 8);
    matrix amp(1, 7);
    dU = L - R;
    dW = Conserved2Primitive(L) - Conserved2Primitive(R);
    dW(1, 2) = (dU(1, 2) - u*dU(1, 1))/rho;
    dW(1, 3) = (dU(1, 3) - v*dU(1, 1))/rho;
    dW(1, 4) = (dU(1, 4) - w*dU(1, 1))/rho;
    dW(1, 5) = (mygamma-1.0)*((0.5*V2-X)*dW(1, 1) - (u*dU(1, 2) + v*dU(1, 3) + w*dU(1, 4)) + dU(1, 5) - (Bx*dU(1, 6) + By*dU(1, 7) + Bz*dU(1, 8))); 
 
    amp(1, 1) = 0.5*(alphaf*(X*dW(1, 1) + dW(1, 5)) + rho*alphas*cs*S*(betay*dW(1, 3) + betaz*dW(1, 4)) - rho*alphaf*cf*dW(1, 2) + rho_sqrt*alphas*a*(betay*dW(1, 7) + betaz*dW(1,8))); 
    amp(1, 2) = 0.5*(betay*dW(1, 4) - betaz*dW(1 ,3) + S*(betay*dW(1, 8)- betaz*dW(1, 7))/rho_sqrt); 
    amp(1, 3) = 0.5*(alphas*(X*dW(1, 1) + dW(1, 5)) - rho*alphaf*cf*S*(betay*dW(1, 3) + betaz*dW(1, 4)) - rho*alphas*cs*dW(1, 2) + rho_sqrt*alphaf*a*(betay*dW(1, 7) + betaz*dW(1,8)));
    amp(1, 4) = (a*a - X)*dW(1, 1) - (P_L - P_R);
    amp(1, 5) = 0.5*(alphas*(X*dW(1, 1) + dW(1, 5)) + rho*alphaf*cf*S*(betay*dW(1, 3) + betaz*dW(1, 4)) + rho*alphas*cs*dW(1, 2) + rho_sqrt*alphaf*a*(betay*dW(1, 7) + betaz*dW(1,8)));
    amp(1, 6) = 0.5*(betaz*dW(1, 3) - betay*dW(1 ,4) + S*(betay*dW(1, 8)- betaz*dW(1, 7))/rho_sqrt);
    amp(1, 7) = 0.5*(alphaf*(X*dW(1, 1) + dW(1, 5)) - rho*alphas*cs*S*(betay*dW(1, 3) + betaz*dW(1, 4)) + rho*alphaf*cf*dW(1, 2) + rho_sqrt*alphas*a*(betay*dW(1, 7) + betaz*dW(1,8)));
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
    EigVec(1, 5) = ( alphaf*(H - B2/rho - u*cf) + alphas*cs*S*(v*betay + w*betaz) - alphas*a*pow(By*By + Bz*Bz, 0.5)/rho_sqrt )/(a*a);
    EigVec(1, 6) = 0.0;
    EigVec(1, 7) = alphas*betay/(rho_sqrt*a);
    EigVec(1, 8) = alphas*betaz/(rho_sqrt*a);
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
    EigVec(3, 5) = ( alphas*(H - B2/rho - u*cs) - alphaf*cf*S*(v*betay + w*betaz) - alphaf*a*pow(By*By + Bz*Bz, 0.5)/rho_sqrt )/(a*a);
    EigVec(3, 6) = 0.0;
    EigVec(3, 7) = -alphaf*betay/(rho_sqrt*a);
    EigVec(3, 8) = -alphaf*betaz/(rho_sqrt*a);
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
    EigVec(5, 5) = ( alphas*(H - B2/rho + u*cs) + alphaf*cf*S*(v*betay + w*betaz) - alphaf*a*pow(By*By + Bz*Bz, 0.5)/rho_sqrt )/(a*a);
    EigVec(5, 6) = 0.0;
    EigVec(5, 7) = -alphaf*betay/(rho_sqrt*a);
    EigVec(5, 8) = -alphaf*betaz/(rho_sqrt*a);
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
    EigVec(7, 5) = ( alphaf*(H - B2/rho + u*cf) - alphas*cs*S*(v*betay + w*betaz) - alphas*a*pow(By*By + Bz*Bz, 0.5)/rho_sqrt )/(a*a);
    EigVec(7, 6) = 0.0;
    EigVec(7, 7) = alphas*betay/(rho_sqrt*a);
    EigVec(7, 8) = alphas*betaz/(rho_sqrt*a);
    
    // compute the fluxes of the left and right states
    matrix flux_L(1, 8);
    matrix flux_R(1, 8);

    flux_L = Conserved2FluxX(L);
    flux_R = Conserved2FluxX(R);

    // compute the Roe flux    
    amp = dotprod(amp, EigVal);
    flux = 0.5 * (flux_L + flux_R) - 0.5 * amp * EigVec;
    return flux;
}

matrix HLLD(const matrix &L, const matrix &R, double Bx)
{
    double rhoL = L(1, 1);
    double rhoR = R(1, 1);
    double uL  = L(1, 2)/L(1, 1);
    double uR  = R(1, 2)/R(1, 1);
    double vL  = L(1, 3)/L(1, 1);
    double vR  = R(1, 3)/R(1, 1);
    double wL  = L(1, 4)/L(1, 1);
    double wR  = R(1, 4)/R(1, 1);
    double pL  = ComputePressure(L);
    double pR  = ComputePressure(R);
    double pTL = ComputePressureStar(L);
    double pTR = ComputePressureStar(R);
    double EL  = L(1, 5);
    double ER  = R(1, 5);
    double ByL = L(1, 7);
    double ByR = R(1, 7);
    double BzL = L(1, 8);
    double BzR = R(1, 8);
     
    double cfL = sqrt((mygamma*pL + (Bx*Bx + ByL*ByL + BzL*BzL) + sqrt(pow((mygamma*pL + (Bx*Bx + ByL*ByL + BzL*BzL)), 2) - 4*mygamma*pL*Bx*Bx))/(2*rhoL));
    double cfR = sqrt((mygamma*pR + (Bx*Bx + ByR*ByR + BzR*BzR) + sqrt(pow((mygamma*pR + (Bx*Bx + ByR*ByR + BzR*BzR)), 2) - 4*mygamma*pR*Bx*Bx))/(2*rhoR));
    double SL = min(uL, uR) - max(cfL, cfR);
    double SR = max(uL, uR) + max(cfL, cfR);

    // compute HLL parameters
    double pTstar = ((SR - uR)*rhoR*pTL - (SL - uL)*rhoL*pTR + rhoL*rhoR*(SR - uR)*(SL - uL)*(uR - uL))/((SR - uR)*rhoR - (SL - uL)*rhoL);
    double SM = (pTR - pTL + rhoL*uL*(SL - uL) - rhoR*uR*(SR - uR))/(rhoL*(SL - uL)-rhoR*(SR - uR));
    double rhostarL = rhoL*(SL - uL)/(SL - SM);
    double rhostarR = rhoR*(SR - uR)/(SR - SM);
    double SstarL = SM - abs(Bx)/sqrt(rhostarL);
    double SstarR = SM + abs(Bx)/sqrt(rhostarR);
    double vstarL = vL - Bx*ByL*(SM - uL)/(rhoL*(SL - uL)*(SL - SM) - Bx*Bx);
    double vstarR = vR - Bx*ByR*(SM - uR)/(rhoR*(SR - uR)*(SR - SM) - Bx*Bx);
    double BystarL = ByL*(rhoL*pow((SL - uL), 2.0) - Bx*Bx)/(rhoL*(SL - uL)*(SL - SM) - Bx*Bx);
    double BystarR = ByR*(rhoR*pow((SR - uR), 2.0) - Bx*Bx)/(rhoR*(SR - uR)*(SR - SM) - Bx*Bx);
    double wstarL = wL - Bx*BzL*(SM - uL)/(rhoL*(SL - uL)*(SL - SM) - Bx*Bx);
    double wstarR = wR - Bx*BzR*(SM - uR)/(rhoR*(SR - uR)*(SR - SM) - Bx*Bx);
    double BzstarL = BzL*(rhoL*pow((SL - uL), 2.0) - Bx*Bx)/(rhoL*(SL - uL)*(SL - SM) - Bx*Bx);
    double BzstarR = BzR*(rhoR*pow((SR - uR), 2.0) - Bx*Bx)/(rhoR*(SR - uR)*(SR - SM) - Bx*Bx);
    double EstarL = ((SL - uL)*EL - pTL*uL + pTstar*SM + Bx*(uL*Bx + vL*ByL + wL*BzL - SM*Bx - vstarL*BystarL - wstarL*BzstarL))/(SL - SM);
    double EstarR = ((SR - uR)*ER - pTR*uR + pTstar*SM + Bx*(uR*Bx + vR*ByR + wR*BzR - SM*Bx - vstarR*BystarR - wstarR*BzstarR))/(SR - SM);
    double rhostar2L = rhostarL;
    double rhostar2R = rhostarR;
    double vstar2 = (sqrt(rhostarL)*vstarL + sqrt(rhostarR)*vstarR + (BystarR - BystarL)*copysign(1.0, Bx))/(sqrt(rhostarL) + sqrt(rhostarR));
    double wstar2 = (sqrt(rhostarL)*wstarL + sqrt(rhostarR)*wstarR + (BzstarR - BzstarL)*copysign(1.0, Bx))/(sqrt(rhostarL) + sqrt(rhostarR));
    double Bystar2 = (sqrt(rhostarL)*BystarR + sqrt(rhostarR)*BystarL + sqrt(rhostarL)*sqrt(rhostarR)*(vstarR-vstarL)*copysign(1.0, Bx))/(sqrt(rhostarL) + sqrt(rhostarR));
    double Bzstar2 = (sqrt(rhostarL)*BzstarR + sqrt(rhostarR)*BzstarL + sqrt(rhostarL)*sqrt(rhostarR)*(wstarR-wstarL)*copysign(1.0, Bx))/(sqrt(rhostarL) + sqrt(rhostarR));
    double Estar2L = EstarL - sqrt(rhostarL)*(SM*Bx + vstarL*BystarL + wstarL*BzstarL - SM*Bx - vstar2*Bystar2 - wstar2*Bzstar2)*copysign(1.0, Bx);
    double Estar2R = EstarR + sqrt(rhostarR)*(SM*Bx + vstarR*BystarR + wstarR*BzstarR - SM*Bx - vstar2*Bystar2 - wstar2*Bzstar2)*copysign(1.0, Bx);
    // compute F
    matrix FL(1, 8);
    matrix FR(1, 8);
    FL = Conserved2FluxX(L);
    FR = Conserved2FluxX(R);
    // compute U*
    matrix UstarL(1, 8);
    matrix UstarR(1, 8);

    UstarL(1, 1) = rhostarL;
    UstarL(1, 2) = rhostarL*SM;
    UstarL(1, 3) = rhostarL*vstarL;
    UstarL(1, 4) = rhostarL*wstarL;
    UstarL(1, 5) = EstarL;
    UstarL(1, 6) = Bx;
    UstarL(1, 7) = BystarL;
    UstarL(1, 8) = BzstarL;

    UstarR(1, 1) = rhostarR;
    UstarR(1, 2) = rhostarR*SM;
    UstarR(1, 3) = rhostarR*vstarR;
    UstarR(1, 4) = rhostarR*wstarR;
    UstarR(1, 5) = EstarR;
    UstarR(1, 6) = Bx;
    UstarR(1, 7) = BystarR;
    UstarR(1, 8) = BzstarR;
    // compute F*
    matrix FstarL(1, 8);
    matrix FstarR(1, 8);
    for(int i = 1; i <= 8; i++)
    {
        FstarL(1, i) = FL(1, i) + SL*(UstarL(1, i) - L(1, i));
        FstarR(1, i) = FR(1, i) + SR*(UstarR(1, i) - R(1, i));
    }
    // compute U**
    matrix Ustar2L(1, 8);
    matrix Ustar2R(1, 8);

    Ustar2L(1, 1) = rhostar2L;
    Ustar2L(1, 2) = rhostar2L*SM;
    Ustar2L(1, 3) = rhostar2L*vstar2;
    Ustar2L(1, 4) = rhostar2L*wstar2;
    Ustar2L(1, 5) = Estar2L;
    Ustar2L(1, 6) = Bx;
    Ustar2L(1, 7) = Bystar2;
    Ustar2L(1, 8) = Bzstar2;
    
    Ustar2R(1, 1) = rhostar2R;
    Ustar2R(1, 2) = rhostar2R*SM;
    Ustar2R(1, 3) = rhostar2R*vstar2;
    Ustar2R(1, 4) = rhostar2R*wstar2;
    Ustar2R(1, 5) = Estar2R;
    Ustar2R(1, 6) = Bx;
    Ustar2R(1, 7) = Bystar2;
    Ustar2R(1, 8) = Bzstar2;
    // compute F**
    matrix Fstar2L(1, 8);
    matrix Fstar2R(1, 8);
    for(int i = 1; i <= 8; i++)
    {
        Fstar2L(1, i) = FstarL(1, i) + SstarL*(Ustar2L(1, i) - UstarL(1, i));
        Fstar2R(1, i) = FstarR(1, i) + SstarR*(Ustar2R(1, i) - UstarR(1, i));
    }

    // compute flux
    matrix F(1, 8);
    if(SL >= 0)
    {
        for(int i = 1; i <= 8; i++)
        {
            F(1, i) = FL(1, i);
        }
    }
    else if(SL < 0 && SstarL >= 0)
    {
        for(int i = 0; i < 8; i++)
        {
             F(1, i) = FstarL(1, i);
        }
    }
    else if(SstarL < 0 && SM >=0)
    {
        for(int i = 0; i < 8; i++)
        {
            F(1, i) = Fstar2L(1, i);
        }
    }
    else if(SM < 0 && SstarR >= 0)
    {
        for(int i = 0; i < 8; i++)
        {
            F(1, i) = Fstar2R(1, i);
        }
    }
    else if(SstarR < 0 && SR >= 0)
    {
        for(int i = 0; i < 8; i++)
        {
            F(1, i) = FstarR(1, i);
        }
    }
    else if(SR < 0)
    {
        for(int i = 0; i < 8; i++)
        {
            F(1, i) = FR(1, i);
        }
    }
    return F;
}

#endif //__SUBFUNC_CPP__
