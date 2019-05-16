#include <iterator>
#include <vector> 
#include <string>
#include <stdio.h>
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <assert.h>
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;
using namespace std;
/*-------------------------------------------------------------------------------
 parameters
-------------------------------------------------------------------------------*/
// constants
const double pi       = 3.14159265358979323846;
const double L        = 1.0;       // 1-D computational domain size
const int N_In        = 64;        // number of computing cells
const double cfl      = 0.8;       // Courant factor
const int nghost      = 2;         // number of ghost zones
const double cs       = 1.0;       // sound speed
const double d_amp    = 1.0e-6;    // density perturbation amplitude
const double d0       = 1.0;       // density background
const double g        = 5.0/3.0;   // ratio of specific heats
const double end_time = 5.0;       // simulation time

// derived constants
const int N  = N_In + 2*nghost;    // total number of cells including ghost zones
const double dx = L/N_In;          // spatial resolution

// plotting parameters
const int nstep_per_image = 1;     // plotting frequency

// global variables
double t = 0.0;
  
/*-------------------------------------------------------------------------------
define the reference solution
-------------------------------------------------------------------------------*/
void ref_func( double x, double *U ){
   double WaveK = 2.0*pi/L;           // wavenumber
   double WaveW = 2.0*pi/(L/cs);      // angular frequency
   double Phase = WaveK*x - WaveW*t;  // wave phase

   double u1 = cs*d_amp/d0;           // velocity perturbation
   double P0 = cs*cs*d0/g;            // background pressure
   double P1 = cs*cs*d_amp;           // pressure perturbation

// d/u/v/w/P/e = density/velocity-x/velocity-y/velocity-z/pressure/total energy
   double d = d0 + d_amp*cos(Phase);
   double u = u1*cos(Phase);
   double v = 0.0;
   double w = 0.0;
   double P = P0 + P1*cos(Phase);
   double E = P/(g-1.0) + 0.5*d*( u*u + v*v + w*w );
   U[0] = d;
   U[1] = d*u;
   U[2] = d*v;
   U[3] = d*w;
   U[4] = E; 
// conserved variables [0/1/2/3/4] <--> [density/momentum x/momentum y/momentum z/energy]
}
/*-------------------------------------------------------------------------------
define boundary condition by setting ghost zones
-------------------------------------------------------------------------------*/
void BoundaryCondition( double U[][5] ){
// periodic
   for(int i = 0; i < nghost; i++){
      for(int j = 0; j < 5; j++){
         U[i][j]   = U[i+N_In][j];
      }
   }
   for(int i = N-nghost; i < N; i++){
      for(int j = 0; j < 5; j++){
         U[i][j] = U[i-N_In][j];
      }
   }
}

/*-------------------------------------------------------------------------------
compute pressure
-------------------------------------------------------------------------------*/
double ComputePressure( double d, double px, double py, double pz, double e ){
   double P = (g-1.0)*( e - 0.5*(px*px + py*py + pz*pz)/d );
   return P;
}

/*-------------------------------------------------------------------------------
compute time-step by the CFL condition
-------------------------------------------------------------------------------*/
double ComputeTimestep( double U[][5] ){
   double P[N];
   double a[N];
   double u[N];
   double v[N];
   double w[N];
   double max_info_speed = 0;
   for(int i = 0; i < N; i++){
      P[i] = ComputePressure( U[i][0], U[i][1], U[i][2], U[i][3], U[i][4] );
      a[i] = pow(( g*P[i]/U[i][0] ), 0.5);
      u[i] = abs( U[i][1]/U[i][0] );
      v[i] = abs( U[i][2]/U[i][0] );
      w[i] = abs( U[i][3]/U[i][0] );
      if( u[i]+a[i] > max_info_speed){
         max_info_speed = u[i]+a[i]; 
      }   
   }
//maximum information speed in 3D
   double dt_cfl = cfl*dx/max_info_speed;
   double dt_end = end_time - t;

   return min( dt_cfl, dt_end );
}

/*-------------------------------------------------------------------------------
compute limited slope
-------------------------------------------------------------------------------*/
double* ComputeLimitedSlope( double *L, double *C, double *R ){
   double slope_L[5];
   double slope_R[5];
   double slope_LR[5];
   static double slope_limited[5];
   for(int i = 0; i < 5; i++){
// compute the left and right slopes
      slope_L[i] = C[i] - L[i];
      slope_R[i] = R[i] - C[i];

// apply the van-Leer limiter
      slope_LR[i]      = slope_L[i]*slope_R[i];
      if(slope_LR[i] > 0.0){
         slope_limited[i] = 2.0*slope_LR[i]/(slope_L[i]+slope_R[i]);    
      }
      else{
         slope_limited[i] = 0;
      }
   }
   return slope_limited;
}

/*-------------------------------------------------------------------------------
convert conserved variables to primitive variables
-------------------------------------------------------------------------------*/
void Conserved2Primitive( double *U, double *W ){
   double A[5]; 
   A[0] = U[0];
   A[1] = U[1]/U[0];
   A[2] = U[2]/U[0];
   A[3] = U[3]/U[0];
   A[4] = ComputePressure( U[0], U[1], U[2], U[3], U[4] );
   for(int i = 0; i < 5; i++){
      W[i] = A[i];
   }
}
/*-------------------------------------------------------------------------------
convert primitive variables to conserved variables
-------------------------------------------------------------------------------*/
void Primitive2Conserved( double *W, double *U){ 
   double A[5]; 
   A[0] = W[0];
   A[1] = W[0]*W[1];
   A[2] = W[0]*W[2];
   A[3] = W[0]*W[3];
   A[4] = W[4]/(g-1.0) + 0.5*W[0]*( W[1]*W[1] + W[2]*W[2] + W[3]*W[3] );
   for(int i = 0; i < 5; i++){
      U[i] = A[i];
   }
}

/*-------------------------------------------------------------------------------
piecewise-linear data reconstruction
-------------------------------------------------------------------------------*/
void DataReconstruction_PLM( double U[][5], double L[][5], double R[][5]){

// allocate memory
   double W[N][5];
   

// conserved variables --> primitive variables
   for (int j = 0; j < N; j++ ){
      Conserved2Primitive( U[j], W[j] );
   }
   for (int j = 1; j < N-1; j++ ){
//    compute the left and right states of each cell
      double *slope_limited = ComputeLimitedSlope( W[j-1], W[j], W[j+1] );

//    get the face-centered variablesi
      for(int i = 0; i < 5; i++){
         L[j][i] = W[j][i] - 0.5*slope_limited[i];
         R[j][i] = W[j][i] + 0.5*slope_limited[i];
      
//    ensure face-centered variables lie between nearby volume-averaged (~cell-centered) values
         L[j][i] = max( L[j][i], min( W[j-1][i], W[j][i] ) );
         L[j][i] = min( L[j][i], max( W[j-1][i], W[j][i] ) );
         R[j][i] = 2.0*W[j][i] - L[j][i];

         R[j][i] = max( R[j][i], min( W[j+1][i], W[j][i] ) );
         R[j][i] = min( R[j][i], max( W[j+1][i], W[j][i] ) );
         L[j][i] = 2.0*W[j][i] - R[j][i];
      }
//    primitive variables --> conserved variables
      Primitive2Conserved( L[j], L[j] );
      Primitive2Conserved( R[j], R[j] );
   }
   
}

/*-------------------------------------------------------------------------------
convert conserved variables to fluxes
-------------------------------------------------------------------------------*/
void Conserved2Flux( double *U, double *flux ){

   double P = ComputePressure( U[0], U[1], U[2], U[3], U[4] );
   double u = U[1] / U[0];

   flux[0] = U[1];
   flux[1] = u*U[1] + P;
   flux[2] = u*U[2];
   flux[3] = u*U[3];
   flux[4] = u*( U[4] + P );
}

/*-------------------------------------------------------------------------------
Roe's Riemann solver
-------------------------------------------------------------------------------*/
void Roe( double *L, double *R, double *flux ){
// compute the enthalpy of the left and right states: H = (E+P)/rho
   double P_L = ComputePressure( L[0], L[1], L[2], L[3], L[4] );
   double P_R = ComputePressure( R[0], R[1], R[2], R[3], R[4] );
   double H_L = ( L[4] + P_L )/L[0];
   double H_R = ( R[4] + P_R )/R[0];

// compute Roe average values
   double rhoL_sqrt = pow(L[0], 0.5);
   double rhoR_sqrt = pow(R[0], 0.5);

   double u  = ( L[1]/rhoL_sqrt + R[1]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt );
   double v  = ( L[2]/rhoL_sqrt + R[2]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt );
   double w  = ( L[3]/rhoL_sqrt + R[3]/rhoR_sqrt ) / ( rhoL_sqrt + rhoR_sqrt );
   double H  = ( rhoL_sqrt*H_L  + rhoR_sqrt*H_R  ) / ( rhoL_sqrt + rhoR_sqrt );
   double V2 = u*u + v*v + w*w;
// check negative pressure
   assert(H-0.5*V2 > 0.0 && "negative pressure!");
   double a = pow( (g-1.0)*(H - 0.5*V2) ,0.5);

// compute the amplitudes of different characteristic waves
   double dU[5];    
   for(int i =0; i < 5; i++){  
      dU[i] = R[i] - L[i];
   }
   double amp[5];
   amp[2] = dU[2] - v*dU[0];
   amp[3] = dU[3] - w*dU[0];
   amp[1] = ((g-1.0)/(a*a))*( dU[0]*(H-u*u) + u*dU[1] - dU[4] + v*amp[2] + w*amp[3] );
   amp[0] = (0.5/a)*( dU[0]*(u+a) - dU[1] - a*amp[1] );
   amp[4] = dU[0] - amp[0] - amp[1];

// compute the eigenvalues and right eigenvector matrix
   double EigenValue[5]    = {u-a, u, u, u, u+a};
   double EigenVector_R[5][5] = { {1.0, u-a,   v,   w,  H-u*a},
                                  {1.0,   u,   v,   w, 0.5*V2},
                                  {0.0, 0.0, 1.0, 0.0,      v},
                                  {0.0, 0.0, 0.0, 1.0,      w},
                                  {1.0, u+a,   v,   w,  H+u*a} };

// compute the fluxes of the left and right states
   double flux_L[5];
   double flux_R[5];
   Conserved2Flux( L, flux_L );
   Conserved2Flux( R, flux_R );
// compute the Roe flux
   for(int i = 0; i < 5; i++){
      amp[i] *= abs( EigenValue[i] );
   }
   for(int i = 0; i < 5; i++){
      double s = 0;
      for(int j = 0; j < 5; j++){
         s = s - 0.5*amp[j]*EigenVector_R[j][i];
      }
      flux[i] = 0.5*( flux_L[i] + flux_R[i] ) + s;
   }
}


/*-------------------------------------------------------------------------------
update animation
-------------------------------------------------------------------------------*/
void update( double *x, double U[][5], double U_ref[][5] ){
   
   for (int i = 0; i < nstep_per_image; i++ ){
      
//    set the boundary conditions
      BoundaryCondition( U );

//    estimate time-step from the CFL condition
      double dt = ComputeTimestep( U );
      printf( "t = %13.7e --> %13.7e, dt = %13.7e  \n", t, t+dt, dt );

//    data reconstruction
      double L[N][5];
      double R[N][5];
      DataReconstruction_PLM( U, L, R );

//    update the face-centered variables by 0.5*dt
      double flux_L[5];
      double flux_R[5];      
      double dflux[5];
      for (int j = 1; j < N-1; j++ ){
         Conserved2Flux( L[j], flux_L);
         Conserved2Flux( R[j], flux_R);
         for(int k = 0; k < 5; k++){
            dflux[k]  = 0.5*(dt/dx)*( flux_R[k] - flux_L[k] );
            L[j][k] -= dflux[k];
            R[j][k] -= dflux[k];
         }
      }
//    compute fluxes
      double flux[N][5];
      for (int j = nghost; j < N-nghost+1; j++ ){
//       R[j-1] is the LEFT state at the j+1/2 inteface
         Roe( R[j-1], L[j], flux[j] );
      }
//    update the volume-averaged input variables by dt
      for(int k = nghost; k < N-nghost; k++){
         for(int j = 0; j < 5; j++){
            U[k][j] -= (dt/dx)*( flux[k+1][j] - flux[k][j] );
         }
      }
//    update time
      t = t + dt;  
   }
   
// calculate the reference analytical solution and estimate errors
   for (int j = 0; j < N_In; j++ ){
      ref_func( x[j], U_ref[j+nghost] );
   }
      
}

/*-------------------------------------------------------------------------------
main
-------------------------------------------------------------------------------*/
int main(){
   double x[N_In];
   double U[N][5];
   double U_ref[N][5];
   double d[N_In];    
   double d_ref[N_In];
   double err = 0;   
 
   for (int j = 0; j < N_In; j++){ 
      x[j] = (j+0.5)*dx;              // cell-centered coordinates
      ref_func( x[j], U[j+nghost] );
      ref_func( x[j], U_ref[j+nghost] );
      d[j] = U[j+nghost][0];
      d_ref[j] = U_ref[j+nghost][0];   
   }
   vector<double> xx(begin(x), end(x));
   vector<double> dd(begin(d), end(d));
   vector<double> dd_ref(begin(d_ref), end(d_ref)); 
   

//   string s = "t = "  t  " err = "  err;
//   plt::title(s);

   plt::xlim( 0.0, L );
   plt::ylim( d0 -1.5*d_amp, d0 +1.5*d_amp );
       
   plt::Plot numerical("numerical", "r-");
   plt::Plot reference("reference", "b--");
   plt::legend();
   numerical.update(xx, dd);
   reference.update(xx, dd_ref);
      
   while(t < end_time){
      update( x, U, U_ref); 
      for (int j = 0; j < N_In; j++ ){
         d[j] = U[j+nghost][0];
         d_ref[j] = U_ref[j+nghost][0];
      }
      dd.clear();
      dd_ref.clear();
      dd.assign(begin(d), end(d));
      dd_ref.assign(begin(d_ref), end(d_ref)); 
//create figure
   
//      plt::plot(x, d, "r--");
//      plt::plot(x, d_ref, "b-");
//      plt::named_plot("numerical", x, d);
//      plt::named_plot("reference", x, d_ref);


//      plt::xlim( 0.0, L );
//      plt::ylim( -1.5*d_amp, +1.5*d_amp );
//      plt::show();

//    plt::title("t = %13.7e err =  %13.7e", t, err);
//      plt::legend();

      numerical.update(xx, dd);
      reference.update(xx, dd_ref);
      plt::pause(0.1);
   }

   return 0;
}



