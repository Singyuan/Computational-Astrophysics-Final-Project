#ifndef __DATARECONSTRUCTION_CPP__
#define __DATARECONSTRUCTION_CPP__

#include "ConstPara.h"

// compute limited slope of one variable
matrix ComputeSlopeSession(matrix d)
{
    double slope_L;
    double slope_R;
    double slope_LR;
    matrix temp(1, N);

    for (int i = 2; i < N; i++)
    {
        slope_L = d(1, i) - d(1, i - 1);
        slope_R = d(1, i+1) - d(1, i);
        slope_LR = slope_L * slope_R;
        if (slope_LR > 0)
        {
            temp(1, i) = 2.0 * slope_LR / (slope_L + slope_R);
        }
        else
        {
            temp(1, i) = 0.0;
        }
    }
    return temp;
}

// compute limited slope
physvar ComputeLimitedSlope(physvar W)
{
    physvar slope_limited;
    slope_limited.d = ComputeSlopeSession(W.d);
    slope_limited.u = ComputeSlopeSession(W.u);
    slope_limited.E = ComputeSlopeSession(W.E);
    return slope_limited;
}

// get the face-centered variables
void FaceCenterVar(matrix d, matrix slope, matrix &tempLd, matrix &tempRd)
{
    for (int i = 2; i < N; i++)
    {
        tempLd(1, i) = d(1, i) - 0.5 * slope(1, i);
        tempRd(1, i) = d(1, i) + 0.5 * slope(1, i);
    }
}

// piecewise-linear data reconstruction
void DataReconstruction_PLM(physvar U, physvar &L, physvar &R)
{
    // allocate memory
    physvar W;
    physvar slope;

    // conserved variables-- > primitive variables
    W = Conserved2Primitive(U);

    // compute the left and right states of each cell
    slope = ComputeLimitedSlope(W);

    //  get the face-centered variables
    matrix tempLd(1, N);
    matrix tempRd(1, N);
    FaceCenterVar(W.d, slope.d, tempLd, tempRd);
    L.d = tempLd;
    R.d = tempRd;

    FaceCenterVar(W.u, slope.u, tempLd, tempRd);
    L.u = tempLd;
    R.u = tempRd;

    FaceCenterVar(W.E, slope.E, tempLd, tempRd);
    L.E = tempLd;
    R.E = tempRd;



    // ensure face-centered variables lie between nearby volume-averaged (~cell-centered) values

    // primitive variables --> conserved variables
    L = Primitive2Conserved(L);
    R = Primitive2Conserved(R);
}

#endif // __DATARECONSTRUCTION_CPP__