#ifndef __CONSTPAPA_H__
#define __CONSTPAPA_H__

#include <iostream>

// #define PI = 3.14159265358

// constants
const double L = 1.0;           //1 - D computational domain size
const int N_In = 10;            //number of computing cells
const double cfl = 0.8;         //Courant factor
const int nghost = 2;           //number of ghost zones
const double gamma = 5.0 / 3.0; //ratio of specific heats
const double end_time = 0.1;    //simulation time
const int largenumber = 1;      //stop the loop

// derived constants
int N = N_In + 2 * nghost; //total number of cells including ghost zones
double dx = L / N_In;      //spatial resolution

#endif // _CONSTPAPA_H__