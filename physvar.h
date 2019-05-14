#ifndef __PHYSVAR_H__
#define __PHYSVAR_H__


#include "matrix.h"
#include "ConstPara.h"

class physvar
{
public:
    physvar(void);
    matrix d; // density
    matrix u; // velocity of x
    matrix E; // energy density
    // int a[5];

private:
    int c;
};

physvar::physvar(void) : d(1, N), u(1, N), E(1, N) {}

#endif // __HYSVAR_H__