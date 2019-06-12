// @File   : matrix.h
// @Author : Singyuan Yeh (yeh76385@gmail.com)
// @Date   : 2019 / 5 / 14

#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <math.h>
class matrix
{
private:
    int Dim0;      //numeber of physics variable
    int Dim1;      //numeber of x grid
    int Dim2;      //numeber of y grid
    int Dim3;
    int Numel;     //numeber of element of matrix
    double *Elem;  //element of matrix
public:
    matrix(void);                               //contruct matrix by default
    matrix(int Dim0s);
    matrix(int Dim0s, int Dim1s);               //contruct matrix by particular numel and default value 0
    matrix(int Dim0s, int Dim1s, int Dim2s);
    matrix(int Dim0s, int Dim1s, int Dim2s, int Dim3s);
    ~matrix(void);

    void display();                      // display matrix
    void showX(int D2, int D3);         // display matrix
    void showY(int D1, int D3);         // display matrix
    void showZ(int D1, int D2);         // display matrix
    matrix GetRow(int D1, int D2, int D3);                // view the value of R row
    void SetRow(int D1, int D2, int D3, const matrix &M); // set the value of R row
    double &operator()(int D0);
    double &operator()(int D0, int D1);
    double &operator()(int D0, int D1, int D2, int D3);
    const matrix &operator=(const matrix &M); // given
    double *GetElem();

    // public:
    friend matrix dotprod(const matrix &matrix1, const matrix &matrix2); // matrix dot product elementwise
    // friend matrix dotpow(const matrix &matrix, const double val);        // matrix dot power elementwise

    friend matrix operator+(const matrix &matrix1, const matrix &matrix2); // matrix add
    friend matrix operator+(const matrix &matrix1, const double val);      // value add (overloading 1)
    friend matrix operator+(const double val, const matrix &matrix1);      // value add (overloading 2)

    friend matrix operator-(const matrix &matrix1, const matrix &matrix2); // matrix minus
    friend matrix operator-(const matrix &matrix1, const double val);      // value minus (overloading 1)
    friend matrix operator-(const double val, const matrix &matrix1);      // value minus (overloading 1)

    friend matrix operator*(const matrix &matrix1, const matrix &matrix2); // matrix multiply LA
    friend matrix operator*(const matrix &matrix1, const double val);      // value multiply (overloading 1)
    friend matrix operator*(const double val, const matrix &matrix1);      // value multiply (overloading 2)

    // friend matrix operator/(const matrix &matrix1, const matrix &matrix2); // matrix dot divide elementwise
    // friend matrix operator/(const matrix &matrix1, const double val);      // value divide
};

// ##########################
// Contructor
// ##########################

// constructor: contruct matrix by default
matrix::matrix(void)
{
    Dim0 = 0;
    Dim1 = 0;
    Dim2 = 0;
    Dim3 = 0;
    Numel = 0;
    Elem = NULL;
}

// constructor: contruct matrix by particular numel and default value 0

matrix::matrix(int Dim0s)
{
    double val = 0.0;
    Dim0 = Dim0s;
    Dim1 = 1;
    Dim2 = 1;
    Dim3 = 1;
    Numel = Dim0 * Dim1 * Dim2 * Dim3;
    Elem = new double[Numel];
    for (int i = 0; i < Numel; i++)
        Elem[i] = val;
}

matrix::matrix(int Dim0s, int Dim1s)
{
    double val = 0.0;
    Dim0 = Dim0s;
    Dim1 = Dim1s;
    Dim2 = 1;
    Dim3 = 1;
    Numel = Dim0 * Dim1 * Dim2 * Dim3;
    Elem = new double[Numel];
    for (int i = 0; i < Numel; i++)
        Elem[i] = val;
}

matrix::matrix(int Dim0s, int Dim1s, int Dim2s)
{
    double val = 0.0;
    Dim0 = Dim0s;
    Dim1 = Dim1s;
    Dim2 = Dim2s;
    Dim3 = 1;
    Numel = Dim0 * Dim1 * Dim2 * Dim3;
    Elem = new double[Numel];
    for (int i = 0; i < Numel; i++)
        Elem[i] = val;
}

matrix::matrix(int Dim0s, int Dim1s, int Dim2s, int Dim3s)
{
    double val = 0.0;
    Dim0 = Dim0s;
    Dim1 = Dim1s;
    Dim2 = Dim2s;
    Dim3 = Dim3s;
    Numel = Dim0 * Dim1 * Dim2 * Dim3;
    Elem = new double[Numel];
    for (int i = 0; i < Numel; i++)
        Elem[i] = val;
}

// destructor
matrix::~matrix(void)
{
    delete[] Elem;
}

// ######################
// Some function
// ######################

// display
void matrix::display(void)
{
    for (int i = 0; i < Dim1; i++)
    {
        for (int j = 0; j < Dim2; j++)
        {
            for (int k = 0; k < Dim3; k++)
            {
                printf("(Dim1, Dim2, Dim3) = (%d, %d, %d)\n", i + 1, j+1, k+1);
                for (int n = 0; n < Dim0; n++)
                {
                    printf("%6.4f ", Elem[k + Dim3 * (j + Dim2 * (i + Dim1 * n))]);
                }
                printf("\n");
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
}

// display
void matrix::showX(int D2, int D3)
{
    D2 = D2 - 1;
    D3 = D3 - 1;
    for (int i = 0; i < Dim0; i++)
    {
        for (int j = 0; j < Dim1; j++)
        {
            printf("%6.4f ", Elem[D3 + Dim3 * (D2 + Dim2 * (j + Dim1 * i))]);
        }
        printf("\n");
    }
    printf("\n");
}

// display
void matrix::showY(int D1, int D3)
{
    D1 = D1 - 1;
    D3 = D3 - 1;
    for (int i = 0; i < Dim0; i++)
    {
        for (int j = 0; j < Dim2; j++)
        {
            printf("%6.4f ", Elem[D3 + Dim3 * (j + Dim2 * (D1 + Dim1 * i))]);
        }
        printf("\n");
    }
    printf("\n");
}

void matrix::showZ(int D1, int D2)
{
    D1 = D1 - 1;
    D2 = D2 - 1;
    for (int i = 0; i < Dim0; i++)
    {
        for (int j = 0; j < Dim1; j++)
        {
            printf("%6.4f ", Elem[j + Dim3 * (D2 + Dim2 * (D1 + Dim1 * i))]);
        }
        printf("\n");
    }
    printf("\n");
}
// view the value of R row
matrix matrix::GetRow(int D1, int D2, int D3)
{
    matrix temp(Dim0);
    D1 = D1 - 1;
    D2 = D2 - 1;
    D3 = D3 - 1;
    for (int i = 0; i < Dim0; i++)
    {
        temp.Elem[i] = Elem[D3 + Dim3 * (D2 + Dim2 * (D1 + Dim1 * i))];
    }
    return temp;
}

// set the value of R row
void matrix::SetRow(int D1, int D2, int D3, const matrix &M)
{
    D1 = D1 - 1;
    D2 = D2 - 1;
    D3 = D3 - 1;
    for (int j = 0; j < Dim0; j++)
    {
        Elem[D3 + Dim3 * (D2 + Dim2 * (D1 + Dim1 * j))] = M.Elem[j];
    }
}

double *matrix::GetElem()
{
    return Elem;
}

// #######################
// operator overriding
// #######################

// operator =
const matrix &matrix::operator=(const matrix &M)
{
    if (this == &M)
    {
        return *this;
    }
    Dim0 = M.Dim0;
    Dim1 = M.Dim1;
    Dim2 = M.Dim2;
    Dim3 = M.Dim3;
    Numel = M.Numel;
    Elem = new double[Numel];
    for (int i = 0; i < Numel; i++)
        Elem[i] = M.Elem[i];
    return *this;
}

// matrix component define ()
double &matrix::operator()(int D0, int D1, int D2, int D3)
{
    D0 = D0 - 1;
    D1 = D1 - 1;
    D2 = D2 - 1;
    D3 = D3 - 1;
    return Elem[D3 + Dim3 * (D2 + Dim2 * (D1 + Dim1 * D0))];
}

double &matrix::operator()(int D0, int D1)
{
    D0 = D0 - 1;
    D1 = D1 - 1;
    return Elem[D1 + Dim1 * D0];
}

double &matrix::operator()(int D0)
{
    return Elem[D0-1];
}

// #######################
// friend function
// #######################

// friend function dotprod matrix dot product elementwise
matrix dotprod(const matrix &matrix1, const matrix &matrix2)
{
    matrix temp(matrix1.Dim0, matrix1.Dim1, matrix1.Dim2, matrix1.Dim3);
    temp = matrix1;
    if (matrix1.Numel == matrix2.Numel)
    {
        for (int i = 0; i < matrix1.Numel; i++)
            temp.Elem[i] *= matrix2.Elem[i];
        return temp;
    }
    else
    {
        printf("Error! Dimension is not match!");
        exit(1);
    }
}

// // friend function dotpow matrix dot power elementwise
// matrix dotpow(const matrix &matrix, const double val)
// {
//     for (int i = 0; i < matrix.Numel; i++)
//         matrix.Elem[i] = pow(matrix.Elem[i], val);
//     return matrix;
// }

// // #######################
// // friend operator overriding
// // #######################

// friend operator + matrix add
matrix operator+(const matrix &matrix1, const matrix &matrix2)
{
    matrix temp(matrix1.Dim0, matrix1.Dim1, matrix1.Dim2, matrix1.Dim3);
    temp = matrix1;
    if (matrix1.Numel == matrix2.Numel)
    {
        for (int i = 0; i < matrix1.Numel; i++)
            temp.Elem[i] += matrix2.Elem[i];
        return temp;
    }
    else
    {
        printf("Error! Dimension is not match!");
        exit(1);
    }
}

//friend operator + value add (overloading 1)
matrix operator+(const matrix &matrix1, const double val)
{
    matrix Temp(matrix1.Dim0, matrix1.Dim1, matrix1.Dim2, matrix1.Dim3);
    Temp = matrix1;
    for (int i = 0; i < matrix1.Numel; i++)
        Temp.Elem[i] += val;
    return Temp;
}

//friend operator + value add (overloading 2)
matrix operator+(const double val, const matrix &matrix)
{
    return matrix + val;
}

// friend operator - matrix minus
matrix operator-(const matrix &matrix1, const matrix &matrix2)
{
    matrix temp(matrix1.Dim0, matrix1.Dim1, matrix1.Dim2, matrix1.Dim3);
    temp = matrix1;
    if (matrix1.Numel == matrix2.Numel)
    {
        for (int i = 0; i < matrix1.Numel; i++)
            temp.Elem[i] -= matrix2.Elem[i];
        return temp;
    }
    else
    {
        printf("Error! Dimension is not match!");
        exit(1);
    }
}

//friend operator - value minus (overloading 1)
matrix operator-(const matrix &matrix1, const double val)
{
    matrix Temp(matrix1.Dim0, matrix1.Dim1, matrix1.Dim2, matrix1.Dim3);
    Temp = matrix1;
    for (int i = 0; i < matrix1.Numel; i++)
        Temp.Elem[i] -= val;
    return Temp;
}

//friend operator - value minus (overloading 2)
matrix operator-(const double val, const matrix &matrix)
{
    return matrix - val;
}

//friend operator * matrix multiply LA
matrix operator*(const matrix &matrix1, const matrix &matrix2)
{
    matrix Temp(matrix1.Dim0, matrix2.Dim1, 1, 1);
    double temp = 0.0;
    int j, k;
    if (matrix1.Dim2 != 1 && matrix1.Dim3 != 1 && matrix2.Dim2 != 1 && matrix2.Dim3 != 1)
    {
        printf("This is not 2 dimension matrix");
        exit(1);
    }
    if (matrix1.Dim1 == matrix2.Dim0)
    {
        for (int i = 0; i < Temp.Dim0; i++)
        {
            for (j = 0; j < Temp.Dim1; j++)
            {
                temp = 0.0;
                for (k = 0; k < matrix1.Dim1; k++)
                {
                    // temp += matrix1.GetElem(i+1, k+1)*matrix2.Elem[k*matrix2.Col+j];
                    temp += matrix1.Elem[i * matrix1.Dim1 + k] * matrix2.Elem[k * matrix2.Dim1 + j];
                }
                Temp.Elem[i * Temp.Dim1 + j] = temp;
            }
        }

        return Temp;
    }
    else
    {
        printf("Error! Dimension is not match!");
        exit(1);
    }
}

//friend operator * value multiply (overloading 1)
matrix operator*(const matrix &matrix1, const double val)
{
    matrix Temp(matrix1.Dim0, matrix1.Dim1, matrix1.Dim2, matrix1.Dim3);
    Temp = matrix1;
    for (int i = 0; i < matrix1.Numel; i++)
        Temp.Elem[i] = Temp.Elem[i] * val;
    return Temp;
}

// friend operator * value multiply (overloading 2)
matrix operator*(const double val, const matrix &matrix1)
{
    return matrix1 * val;
}

// // friend operator / matrix dot divide elementwise
// matrix operator/(const matrix &matrix1, const matrix &matrix2)
// {
//     matrix temp(matrix1.Row, matrix1.Col);
//     temp = matrix1;
//     if (matrix1.Col == matrix2.Col && matrix1.Row == matrix2.Row)
//     {
//         for (int i = 0; i < matrix1.Numel; i++)
//             temp.Elem[i] /= matrix2.Elem[i];
//         return temp;
//     }
//     else
//     {
//         printf("Error! Dimension is not match!");
//         exit(1);
//     }
// }

// // friend operator / value divide
// matrix operator/(const matrix &matrix1, const double val)
// {
//     for (int i = 0; i < matrix1.Numel; i++)
//         matrix1.Elem[i] = matrix1.Elem[i] / val;
//     return matrix1;
// }

#endif // __MATRIX_H__