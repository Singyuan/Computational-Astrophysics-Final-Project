#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <cmath>
#include <algorithm>

double SpeedX(double rhoL, double rhoR, double uL, double uR, double pL, double pR, double Bx, double ByL, double ByR, double BzL, double BzR, double gamma){
	double pTL=pL+(pow(Bx,2)+pow(ByL,2)+pow(BzL,2))/2;
	double pTR=pR+(pow(Bx,2)+pow(ByR,2)+pow(BzR,2))/2;
	double cfL=sqrt((gamma*pL+(pow(Bx,2)+pow(ByL,2)+pow(BzL,2))+sqrt(pow(gamma*pL+(pow(Bx,2)+pow(ByL,2)+pow(BzL,2)),2)-4*gamma*pL*pow(Bx,2)))/(2*rhoL));
	double cfR=sqrt((gamma*pR+(pow(Bx,2)+pow(ByR,2)+pow(BzR,2))+sqrt(pow(gamma*pR+(pow(Bx,2)+pow(ByR,2)+pow(BzR,2)),2)-4*gamma*pR*pow(Bx,2)))/(2*rhoR));
	double SL=std::min(uL,uR)-std::max(cfL,cfR);
	double SR=std::max(uL,uR)+std::max(cfL,cfR);
	double Sx=std::max(std::abs(SL),std::abs(SR));
	return Sx;
}

double* FluxX(double rhoL, double rhoR, double uL, double uR, double vL, double vR, double wL, double wR, double pL, double pR, double Bx, double ByL, double ByR, double BzL, double BzR, double gamma){
	double pTL=pL+(pow(Bx,2)+pow(ByL,2)+pow(BzL,2))/2;
	double pTR=pR+(pow(Bx,2)+pow(ByR,2)+pow(BzR,2))/2;
	double eL=pL/(gamma-1);
	double EL=rhoL*(pow(uL,2)/2+pow(vL,2)/2+pow(wL,2)/2)+eL+(pow(Bx,2)+pow(ByL,2)+pow(BzL,2))/2;
	double eR=pR/(gamma-1);
	double ER=rhoR*(pow(uR,2)/2+pow(vR,2)/2+pow(wR,2)/2)+eR+(pow(Bx,2)+pow(ByR,2)+pow(BzR,2))/2;
	double cfL=sqrt((gamma*pL+(pow(Bx,2)+pow(ByL,2)+pow(BzL,2))+sqrt(pow((gamma*pL+(pow(Bx,2)+pow(ByL,2)+pow(BzL,2))),2)-4*gamma*pL*pow(Bx,2)))/(2*rhoL));
	double cfR=sqrt((gamma*pR+(pow(Bx,2)+pow(ByR,2)+pow(BzR,2))+sqrt(pow((gamma*pR+(pow(Bx,2)+pow(ByR,2)+pow(BzR,2))),2)-4*gamma*pR*pow(Bx,2)))/(2*rhoR));
	double SL=std::min(uL,uR)-std::max(cfL,cfR);
	double SR=std::max(uL,uR)+std::max(cfL,cfR);
	double pTstar=((SR-uR)*rhoR*pTL-(SL-uL)*rhoL*pTR+rhoL*rhoR*(SR-uR)*(SL-uL)*(uR-uL))/((SR-uR)*rhoR-(SL-uL)*rhoL);
    double SM=(pTR-pTL+rhoL*uL*(SL-uL)-rhoR*uR*(SR-uR))/(rhoL*(SL-uL)-rhoR*(SR-uR));
    double rhostarL=rhoL*(SL-uL)/(SL-SM);
    double rhostarR=rhoR*(SR-uR)/(SR-SM);
    double SstarL=SM-std::abs(Bx)/sqrt(rhostarL);
    double SstarR=SM+std::abs(Bx)/sqrt(rhostarR);
    double vstarL=vL-Bx*ByL*(SM-uL)/(rhoL*(SL-uL)*(SL-SM)-pow(Bx,2));
    double vstarR=vR-Bx*ByR*(SM-uR)/(rhoR*(SR-uR)*(SR-SM)-pow(Bx,2));
    double BystarL=ByL*(rhoL*pow((SL-uL),2)-pow(Bx,2))/(rhoL*(SL-uL)*(SL-SM)-pow(Bx,2));
    double BystarR=ByR*(rhoR*pow((SR-uR),2)-pow(Bx,2))/(rhoR*(SR-uR)*(SR-SM)-pow(Bx,2));
	double wstarL=wL-Bx*BzL*(SM-uL)/(rhoL*(SL-uL)*(SL-SM)-pow(Bx,2));
	double wstarR=wR-Bx*BzR*(SM-uR)/(rhoR*(SR-uR)*(SR-SM)-pow(Bx,2));
	double BzstarL=BzL*(rhoL*pow((SL-uL),2)-pow(Bx,2))/(rhoL*(SL-uL)*(SL-SM)-pow(Bx,2));
    double BzstarR=BzR*(rhoR*pow((SR-uR),2)-pow(Bx,2))/(rhoR*(SR-uR)*(SR-SM)-pow(Bx,2));
    double EstarL=((SL-uL)*EL-pTL*uL+pTstar*SM+Bx*(uL*Bx+vL*ByL+wL*BzL-SM*Bx-vstarL*BystarL-wstarL*BzstarL))/(SL-SM);
    double EstarR=((SR-uR)*ER-pTR*uR+pTstar*SM+Bx*(uR*Bx+vR*ByR+wR*BzR-SM*Bx-vstarR*BystarR-wstarR*BzstarR))/(SR-SM);
    double rhostar2L=rhostarL;
	double rhostar2R=rhostarR;
    double vstar2=(sqrt(rhostarL)*vstarL+sqrt(rhostarR)*vstarR+(BystarR-BystarL)*copysign(1,Bx))/(sqrt(rhostarL)+sqrt(rhostarR));
    double wstar2=(sqrt(rhostarL)*wstarL+sqrt(rhostarR)*wstarR+(BzstarR-BzstarL)*copysign(1,Bx))/(sqrt(rhostarL)+sqrt(rhostarR));
    double Bystar2=(sqrt(rhostarL)*BystarR+sqrt(rhostarR)*BystarL+sqrt(rhostarL)*sqrt(rhostarR)*(vstarR-vstarL)*copysign(1,Bx))/(sqrt(rhostarL)+sqrt(rhostarR));
	double Bzstar2=(sqrt(rhostarL)*BzstarR+sqrt(rhostarR)*BzstarL+sqrt(rhostarL)*sqrt(rhostarR)*(wstarR-wstarL)*copysign(1,Bx))/(sqrt(rhostarL)+sqrt(rhostarR));
	double Estar2L=EstarL-sqrt(rhostarL)*(SM*Bx+vstarL*BystarL+wstarL*BzstarL-SM*Bx-vstar2*Bystar2-wstar2*Bzstar2)*copysign(1,Bx);
	double Estar2R=EstarR+sqrt(rhostarR)*(SM*Bx+vstarR*BystarR+wstarR*BzstarR-SM*Bx-vstar2*Bystar2-wstar2*Bzstar2)*copysign(1,Bx);
	double UL[8]={rhoL, rhoL*uL, rhoL*vL, rhoL*wL, Bx, ByL, BzL, EL};
    double UR[8]={rhoR, rhoR*uR, rhoR*vR, rhoR*wR, Bx, ByR, BzR, ER};
    double FL[8]={rhoL*uL, (rhoL*pow(uL,2)+pTL-pow(Bx,2)), (rhoL*uL*vL-Bx*ByL), (rhoL*uL*wL-Bx*BzL), 0, (ByL*uL-Bx*vL), (BzL*uL-Bx*wL), (uL*(EL+pTL)-Bx*(uL*Bx+vL*ByL+wL*BzL))};
	double FR[8]={rhoR*uR, (rhoR*pow(uR,2)+pTR-pow(Bx,2)), (rhoR*uR*vR-Bx*ByR), (rhoR*uR*wR-Bx*BzR), 0, (ByR*uR-Bx*vR), (BzR*uR-Bx*wR), (uR*(ER+pTR)-Bx*(uR*Bx+vR*ByR+wR*BzR))};
    double UstarL[8]={rhostarL, rhostarL*SM, rhostarL*vstarL, rhostarL*wstarL, Bx, BystarL, BzstarL, EstarL};
    double UstarR[8]={rhostarR, rhostarR*SM, rhostarR*vstarR, rhostarR*wstarR, Bx, BystarR, BzstarR, EstarR};
    double FstarL[8], FstarR[8];
	for(int u=0; u<8; u++){
    	FstarL[u]=FL[u]+SL*(UstarL[u]-UL[u]);
    	FstarR[u]=FR[u]+SR*(UstarR[u]-UR[u]);
	}
    double Ustar2L[8]={rhostar2L, rhostar2L*SM, rhostar2L*vstar2, rhostar2L*wstar2, Bx, Bystar2, Bzstar2, Estar2L};
    double Ustar2R[8]={rhostar2R, rhostar2R*SM, rhostar2R*vstar2, rhostar2R*wstar2, Bx, Bystar2, Bzstar2, Estar2R};
    double Fstar2L[8], Fstar2R[8];
	for(int u=0; u<8; u++){
    	Fstar2L[u]=FstarL[u]+SstarL*(Ustar2L[u]-UstarL[u]);
    	Fstar2R[u]=FstarR[u]+SstarR*(Ustar2R[u]-UstarR[u]);
	}
    double F[8];
    if(SL>=0){
    	for(int u=0; u<8; u++){
    		F[u]=FL[u];
		}
	}
	else if(SL<0 && SstarL>=0){
		for(int u=0; u<8; u++){
    		F[u]=FstarL[u];
		}
	}
	else if(SstarL<0 && SM>=0){
		for(int u=0; u<8; u++){
    		F[u]=Fstar2L[u];
		}
	}
	else if(SM<0 && SstarR>=0){
		for(int u=0; u<8; u++){
    		F[u]=Fstar2R[u];
		}
	}
	else if(SstarR<0 && SR>=0){
		for(int u=0; u<8; u++){
			F[u]=FstarR[u];
		}
	}
	else if(SR<0){
		for(int u=0; u<8; u++){
			F[u]=FR[u];
		}
	}
	return F;
}

double ComputeLimitedSlope(double L, double C, double R){
	// compute the left and right slopes
	double slopeL = C - L;
	double slopeR = R - C;
	double slopelimited;
	// apply the van-Leer limiter
	double slopeLR = slopeL*slopeR;
	if(slopeLR>0.0){
		slopelimited=2.0*slopeLR/(slopeL+slopeR);
	}
	else if(slopeLR<=0.0){
		slopelimited=0.0;
	}
	
	return slopelimited;
}

double reconstructionL(double um, double u, double up){
	double du=ComputeLimitedSlope(um, u, up);
	double L=u+0.5*du;
	L=std::max(L,std::min(u,up));
	L=std::min(L,std::max(u,up));
	return L;
}

double reconstructionR(double um, double u, double up){
	double du=ComputeLimitedSlope(um, u, up);
	double R=u-0.5*du;
	R=std::max(R,std::min(um,u));
	R=std::min(R,std::max(um,u));
	return R;
}

int main(){
	double pi=3.14159265359;
	int Nx=104, Ny=5, Nz=5;
	double dx=1.0/(Nx-4.0), dy=1/(Ny-4.0), dz=1/(Nz-4.0);
	double W[5][Nx][Ny][Nz], E[Nx][Ny][Nz], Bx[Nx][Ny][Nz], By[Nx][Ny][Nz], Bz[Nx][Ny][Nz];
	double Wmarch[5][Nx][Ny][Nz], Emarch[Nx][Ny][Nz], Bxmarch[Nx][Ny][Nz], Bymarch[Nx][Ny][Nz], Bzmarch[Nx][Ny][Nz];
	
	// Initial Condition 
	/*double rhoL=1.08, uL=1.2, vL=0.01, wL=0.5, pL=0.95, ByL=3.6/sqrt(4.0*pi), BzL=2.0/sqrt(4.0*pi);
	double rhoR=1, uR=0, vR=0, wR=0, pR=1, ByR=4/sqrt(4*pi), BzR=2.0/sqrt(4.0*pi);
	int shift=52;
	double tf=0.2, gamma=5.0/3.0, Bxx=4/sqrt(4*pi);*/
	
	double rhoL=1.0, uL=0, vL=0, wL=0, pL=1.0, ByL=1.0, BzL=0;
	double rhoR=0.125, uR=0, vR=0, wR=0, pR=0.1, ByR=-1.0, BzR=0;
	int shift=52;
	double tf=0.1, gamma=5.0/3.0, Bxx=0.75;
	
	for(int m=0; m<Nx; m++){
		double x=(m-shift)*dx;
		for(int n=0; n<Ny; n++){
			for(int l=0; l<Nz; l++){
				if (x<0){
					W[0][m][n][l]=rhoL;
					W[1][m][n][l]=uL;
					W[2][m][n][l]=vL;
					W[3][m][n][l]=wL;
					W[4][m][n][l]=pL;
					Bx[m][n][l]=Bxx;
					By[m][n][l]=ByL;
					Bz[m][n][l]=BzL;
				}
				else{
					W[0][m][n][l]=rhoR;
					W[1][m][n][l]=uR;
					W[2][m][n][l]=vR;
					W[3][m][n][l]=wR;
					W[4][m][n][l]=pR;
					Bx[m][n][l]=Bxx;
					By[m][n][l]=ByR;
					Bz[m][n][l]=BzR;
				}
				E[m][n][l]=W[0][m][n][l]*((pow(W[1][m][n][l],2)+pow(W[2][m][n][l],2)+pow(W[3][m][n][l],2))/2+W[4][m][n][l]/(gamma-1))+(pow(Bx[m][n][l],2)+pow(By[m][n][l],2)+pow(Bz[m][n][l],2))/2;
			}
		}
	}
	
	// START OF MAIN INTEGRATION LOOP
	double t=0, Smax=0;
	while(t<tf){
	//for(int k=0; k<1; k++){
		
		// Compute Time Step
		for(int m=2; m<Nx-1; m++){
			for(int n=2; n<Ny-1; n++){
				for(int l=2; l<Nz-1; l++){
					double Sx=SpeedX(W[0][m-1][n][l], W[0][m][n][l], W[1][m-1][n][l], W[1][m][n][l], W[4][m-1][n][l], W[4][m][n][l], Bx[m][n][l], (By[m-1][n][l]+By[m-1][n+1][l])/2, (By[m][n][l]+By[m][n+1][l])/2, (Bz[m-1][n][l]+Bz[m-1][n][l+1])/2, (Bz[m][n][l]+Bz[m][n][l+1])/2, gamma);
					if (Sx>Smax){
						Smax=Sx;
					}
				}
			}
		}
		double cfl=0.4, d=std::min(dx,dy);
		d=std::min(d,dz);
		double dt=cfl*d/Smax;
		
		// Time Marching
		for(int m=2; m<Nx-2; m++){
			for(int n=2; n<Ny-2; n++){
				for(int l=2; l<Nz-2; l++){
					double U[8]={W[0][m][n][l], W[0][m][n][l]*W[1][m][n][l], W[0][m][n][l]*W[2][m][n][l], W[0][m][n][l]*W[3][m][n][l], Bx[m][n][l], By[m][n][l], Bz[m][n][l], E[m][n][l]};
					double rhoL=reconstructionL(W[0][m-2][n][l], W[0][m-1][n][l], W[0][m][n][l]);
					double rhoR=reconstructionR(W[0][m-1][n][l], W[0][m][n][l], W[0][m+1][n][l]);
					double uL=reconstructionL(W[1][m-2][n][l], W[1][m-1][n][l], W[1][m][n][l]);
					double uR=reconstructionR(W[1][m-1][n][l], W[1][m][n][l], W[1][m+1][n][l]);
					double vL=reconstructionL(W[2][m-2][n][l], W[2][m-1][n][l], W[2][m][n][l]);
					double vR=reconstructionR(W[2][m-1][n][l], W[2][m][n][l], W[2][m+1][n][l]);
					double wL=reconstructionL(W[3][m-2][n][l], W[3][m-1][n][l], W[3][m][n][l]);
					double wR=reconstructionR(W[3][m-1][n][l], W[3][m][n][l], W[3][m+1][n][l]);
					double pL=reconstructionL(W[4][m-2][n][l], W[4][m-1][n][l], W[4][m][n][l]);
					double pR=reconstructionR(W[4][m-1][n][l], W[4][m][n][l], W[4][m+1][n][l]);
					double ByL=reconstructionL(By[m-2][n][l], By[m-1][n][l], By[m][n][l]);
					double ByR=reconstructionR(By[m-1][n][l], By[m][n][l], By[m+1][n][l]);
					double BzL=reconstructionL(Bz[m-2][n][l], Bz[m-1][n][l], Bz[m][n][l]);
					double BzR=reconstructionR(Bz[m-1][n][l], Bz[m][n][l], Bz[m+1][n][l]);
					double *FL=FluxX(rhoL, rhoR, uL, uR, vL, vR, wL, wR, pL, pR, Bx[m][n][l], ByL, ByR, BzL, BzR, gamma);
					for(int u=0; u<8; u++){
						U[u]=U[u]+dt/dx*FL[u];
					}
					rhoL=reconstructionL(W[0][m-1][n][l], W[0][m][n][l], W[0][m+1][n][l]);
					rhoR=reconstructionR(W[0][m][n][l], W[0][m+1][n][l], W[0][m+2][n][l]);
					uL=reconstructionL(W[1][m-1][n][l], W[1][m][n][l], W[1][m+1][n][l]);
					uR=reconstructionR(W[1][m][n][l], W[1][m+1][n][l], W[1][m+2][n][l]);
					vL=reconstructionL(W[2][m-1][n][l], W[2][m][n][l], W[2][m+1][n][l]);
					vR=reconstructionR(W[2][m][n][l], W[2][m+1][n][l], W[2][m+2][n][l]);
					wL=reconstructionL(W[3][m-1][n][l], W[3][m][n][l], W[3][m+1][n][l]);
					wR=reconstructionR(W[3][m][n][l], W[3][m+1][n][l], W[3][m+2][n][l]);
					pL=reconstructionL(W[4][m-1][n][l], W[4][m][n][l], W[4][m+1][n][l]);
					pR=reconstructionR(W[4][m][n][l], W[4][m+1][n][l], W[4][m+2][n][l]);
					ByL=reconstructionL(By[m-1][n][l], By[m][n][l], By[m+1][n][l]);
					ByR=reconstructionR(By[m][n][l], By[m+1][n][l], By[m+2][n][l]);
					BzL=reconstructionL(Bz[m-1][n][l], Bz[m][n][l], Bz[m+1][n][l]);
					BzR=reconstructionR(Bz[m][n][l], Bz[m+1][n][l], Bz[m+2][n][l]);
					double *FR=FluxX(rhoL, rhoR, uL, uR, vL, vR, wL, wR, pL, pR, Bx[m][n][l], ByL, ByR, BzL, BzR, gamma);
					for(int u=0; u<8; u++){
						U[u]=U[u]-dt/dx*FR[u];
					}
					Wmarch[0][m][n][l]=U[0];
					Wmarch[1][m][n][l]=U[1]/U[0];
					Wmarch[2][m][n][l]=U[2]/U[0];
					Wmarch[3][m][n][l]=U[3]/U[0];
					Bxmarch[m][n][l]=U[4];
	                Bymarch[m][n][l]=U[5];
	                Bzmarch[m][n][l]=U[6];
					Wmarch[4][m][n][l]=(U[7]-(pow(Wmarch[1][m][n][l],2)/2+pow(Wmarch[2][m][n][l],2)/2+pow(Wmarch[3][m][n][l],2)/2)*U[0]-(pow(Bxmarch[m][n][l],2)+pow(Bymarch[m][n][l],2)+pow(Bzmarch[m][n][l],2))/2)*(gamma-1.0);
					Emarch[m][n][l]=U[7];
				}
			}
		}
		
		// Boundary Condition
		for(int n=0; n<Ny; n++){
			for(int l=0; l<Nz; l++){
				for(int u=0; u<5; u++){
					Wmarch[u][0][n][l]=Wmarch[u][2][n][l];
					Wmarch[u][Nx-1][n][l]=Wmarch[u][Nx-3][n][l];
				}
				Bxmarch[0][n][l]=Bxmarch[2][n][l];
				Bxmarch[Nx-1][n][l]=Bxmarch[Nx-3][n][l];
				Bymarch[0][n][l]=Bymarch[2][n][l];
				Bymarch[Nx-1][n][l]=Bymarch[Nx-3][n][l];
				Bzmarch[0][n][l]=Bzmarch[2][n][l];
				Bzmarch[Nx-1][n][l]=Bzmarch[Nx-3][n][l];
				Emarch[0][n][l]=Emarch[2][n][l];
				Emarch[Nx-1][n][l]=Emarch[Nx-3][n][l];
				for(int u=0; u<5; u++){
					Wmarch[u][1][n][l]=Wmarch[u][2][n][l];
					Wmarch[u][Nx-2][n][l]=Wmarch[u][Nx-3][n][l];
				}
				Bxmarch[1][n][l]=Bxmarch[2][n][l];
				Bxmarch[Nx-2][n][l]=Bxmarch[Nx-3][n][l];
				Bymarch[1][n][l]=Bymarch[2][n][l];
				Bymarch[Nx-2][n][l]=Bymarch[Nx-3][n][l];
				Bzmarch[1][n][l]=Bzmarch[2][n][l];
				Bzmarch[Nx-2][n][l]=Bzmarch[Nx-3][n][l];
				Emarch[1][n][l]=Emarch[2][n][l];
				Emarch[Nx-2][n][l]=Emarch[Nx-3][n][l];
			}
		}
		for(int m=0; m<Nx; m++){
			for(int l=0; l<Nz; l++){
				for(int u=0; u<5; u++){
					Wmarch[u][m][0][l]=Wmarch[u][m][2][l];
					Wmarch[u][m][Ny-1][l]=Wmarch[u][m][Ny-3][l];
				}
				Bxmarch[m][0][l]=Bxmarch[m][2][l];
				Bxmarch[m][Ny-1][l]=Bxmarch[m][Ny-3][l];
				Bymarch[m][0][l]=Bymarch[m][2][l];
				Bymarch[m][Ny-1][l]=Bymarch[m][Ny-3][l];
				Bzmarch[m][0][l]=Bzmarch[m][2][l];
				Bzmarch[m][Ny-1][l]=Bzmarch[m][Ny-3][l];
				Emarch[m][0][l]=Emarch[m][2][l];
				Emarch[m][Ny-1][l]=Emarch[m][Ny-3][l];
				for(int u=0; u<5; u++){
					Wmarch[u][m][1][l]=Wmarch[u][m][2][l];
					Wmarch[u][m][Ny-2][l]=Wmarch[u][m][Ny-3][l];
				}
				Bxmarch[m][1][l]=Bxmarch[m][2][l];
				Bxmarch[m][Ny-2][l]=Bxmarch[m][Ny-3][l];
				Bymarch[m][1][l]=Bymarch[m][2][l];
				Bymarch[m][Ny-2][l]=Bymarch[m][Ny-3][l];
				Bzmarch[m][1][l]=Bzmarch[m][2][l];
				Bzmarch[m][Ny-2][l]=Bzmarch[m][Ny-3][l];
				Emarch[m][1][l]=Emarch[m][2][l];
				Emarch[m][Ny-2][l]=Emarch[m][Ny-3][l];
			}
		}
		for(int m=0; m<Nx; m++){
			for(int n=0; n<Ny; n++){
				for(int u=0; u<5; u++){
					Wmarch[u][m][n][0]=Wmarch[u][m][n][2];
					Wmarch[u][m][n][Nz-1]=Wmarch[u][m][n][Nz-3];
				}
				Bxmarch[m][n][0]=Bxmarch[m][n][2];
				Bxmarch[m][n][Nz-1]=Bxmarch[m][n][Nz-3];
				Bymarch[m][n][0]=Bymarch[m][n][2];
				Bymarch[m][n][Nz-1]=Bymarch[m][n][Nz-3];
				Bzmarch[m][n][0]=Bzmarch[m][n][2];
				Bzmarch[m][n][Nz-1]=Bzmarch[m][n][Nz-3];
				Emarch[m][n][0]=Emarch[m][n][2];
				Emarch[m][n][Nz-1]=Emarch[m][n][Nz-3];
				for(int u=0; u<5; u++){
					Wmarch[u][m][n][1]=Wmarch[u][m][n][2];
					Wmarch[u][m][n][Nz-2]=Wmarch[u][m][n][Nz-3];
				}
				Bxmarch[m][n][1]=Bxmarch[m][n][2];
				Bxmarch[m][n][Nz-2]=Bxmarch[m][n][Nz-3];
				Bymarch[m][n][1]=Bymarch[m][n][2];
				Bymarch[m][n][Nz-2]=Bymarch[m][n][Nz-3];
				Bzmarch[m][n][1]=Bzmarch[m][n][2];
				Bzmarch[m][n][Nz-2]=Bzmarch[m][n][Nz-3];
				Emarch[m][n][1]=Emarch[m][n][2];
				Emarch[m][n][Nz-2]=Emarch[m][n][Nz-3];
			}
		}
		for(int m=0; m<Nx; m++){
			for(int n=0; n<Ny; n++){
				for(int l=0; l<Nz; l++){
					for(int u=0; u<5; u++){
						W[u][m][n][l]=Wmarch[u][m][n][l];
					}
					Bx[m][n][l]=Bxmarch[m][n][l];
					By[m][n][l]=Bymarch[m][n][l];
					Bz[m][n][l]=Bzmarch[m][n][l];
					E[m][n][l]=Emarch[m][n][l];
				}
			}
		} 
			
		t=t+dt;
		printf("t = %lf, dt = %lf, Smax = %lf, d = %lf\n",t,dt,Smax,d);
	}
	
	// Print The Results 
	for(int m=0; m<Nx; m++){
		printf("%lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf\n",W[0][m][2][2], W[1][m][2][2], W[2][m][2][2], W[3][m][2][2], W[4][m][2][2], Bx[m][2][2], By[m][2][2], Bz[m][2][2], E[m][2][2]);
	}
	
	return 0;
}
