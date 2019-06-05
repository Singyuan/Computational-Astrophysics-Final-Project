#ifndef __CONSTPAPA_H__
#define __CONSTPAPA_H__

// #define PI = 3.14159265358

// constants
const double L = 1.0;             //1-D computational domain size
const int N_In = 64;                 //number of computing cells
const double cfl = 0.3;           //Courant factor
const int nghost = 2;             //number of ghost zones
const double mygamma = 2.0; //ratio of specific heats
const double end_time = 0.1;      //simulation time
const int largenumber = 35;       //stop the loop

// derived constants
int N = N_In + 2 * nghost;    // total number of cells including ghost zones
double dx = L / N_In;         // resolution
const int NThread = 2;        // open mp threads 

#endif // _CONSTPAPA_H__
