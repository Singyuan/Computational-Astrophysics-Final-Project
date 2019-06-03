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
    int Col;      //numeber of column of matrix
    int Row;      //numeber of row of matrix
    int Numel;    //numeber of element of matrix
    double *Elem; //element of matrix
public:
    matrix(void);                               //contruct matrix by default
    matrix(int Rows, int Cols, double Val = 0); //contruct matrix by particular numel and default value 0
    ~matrix(void);

    void display();                      // display matrix
    matrix GetRow(int R);                // view the value of R row
    void SetRow(int R, const matrix &M); // set the value of R row
    matrix GetCol(int C);                // view the value of C col
    void SetCol(int C, const matrix &M); // set the value of C col
    matrix T();                          // transpose

    double &operator()(int R, int C);
    const matrix &operator=(const matrix &M); // given

public:
    friend matrix dotprod(const matrix &matrix1, const matrix &matrix2); // matrix dot product elementwise
    friend matrix dotpow(const matrix &matrix, const double val);        // matrix dot power elementwise

    friend matrix operator+(const matrix &matrix1, const matrix &matrix2); // matrix add
    friend matrix operator+(const matrix &matrix1, const double val);      // value add (overloading 1)
    friend matrix operator+(const double val, const matrix &matrix1);      // value add (overloading 2)

    friend matrix operator-(const matrix &matrix1, const matrix &matrix2); // matrix minus
    friend matrix operator-(const matrix &matrix1, const double val);      // value minus (overloading 1)
    friend matrix operator-(const double val, const matrix &matrix1);      // value minus (overloading 1)

    friend matrix operator*(const matrix &matrix1, const matrix &matrix2); // matrix multiply LA
    friend matrix operator*(const matrix &matrix1, const double val);      // value multiply (overloading 1)
    friend matrix operator*(const double val, const matrix &matrix1);      // value multiply (overloading 2)

    friend matrix operator/(const matrix &matrix1, const matrix &matrix2); // matrix dot divide elementwise
    friend matrix operator/(const matrix &matrix1, const double val);      // value divide
};

// ##########################
// Contructor
// ##########################

// constructor: contruct matrix by default
matrix::matrix(void)
{
    Row = 0;
    Col = 0;
    Numel = Row * Col;
    Elem = NULL;
}

// constructor: contruct matrix by particular numel and default value 0
matrix::matrix(int Rows, int Cols, double val)
{
    Row = Rows;
    Col = Cols;
    Numel = Row * Col;
    Elem = new double[Row * Col];
    for (int i = 0; i < Row * Col; i++)
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

// transpose
matrix matrix::T()
{
    matrix Temp(Col, Row);
    for (int i = 0; i < Temp.Row; i++)
        for (int j = 0; j < Temp.Col; j++)
            Temp.Elem[i * Temp.Col + j] = Elem[j * Temp.Row + i];
    return Temp;
}

// display
void matrix::display(void)
{
    for (int i = 0; i < Row; i++)
    {
        for (int j = 0; j < Col; j++)
        {
            printf("%f ", Elem[i * Col + j]);
        }
        printf("\n");
    }
    printf("\n");
}

// view the value of R row
matrix matrix::GetRow(int R)
{
    matrix temp(1, Col);
#pragma omp parallel for
    for (int j = 0; j < Col; j++)
    {
        temp.Elem[j] = Elem[(R - 1) * Col + j];
    }
    return temp;
}

// set the value of R row
void matrix::SetRow(int R, const matrix &M)
{
#pragma omp parallel for
    for (int j = 0; j < Col; j++)
    {
        Elem[(R - 1) * Col + j] = M.Elem[j];
    }
}

// view the value of C col
matrix matrix::GetCol(int C)
{
    matrix temp(1, Row);
    for (int i = 0; i < Row; i++)
    {
        temp.Elem[i] = Elem[i * Col + (C - 1)];
    }
    return temp;
}

// set the value of C col
void matrix::SetCol(int C, const matrix &M)
{
    for (int i = 0; i < Row; i++)
    {
        Elem[i * Col + (C - 1)] = M.Elem[i];
    }
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
    Col = M.Col;
    Row = M.Row;
    Numel = M.Numel;
    Elem = new double[Numel];

#pragma omp parallel for
    for (int i = 0; i < Numel; i++)
        Elem[i] = M.Elem[i];
    return *this;
}

// matrix component define ()
double &matrix::operator()(int R, int C)
{
    return Elem[(R - 1) * Col + (C - 1)];
}

// #######################
// friend function
// #######################

// friend function dotprod matrix dot product elementwise
matrix dotprod(const matrix &matrix1, const matrix &matrix2)
{
    matrix temp(matrix1.Row, matrix1.Col);
    temp = matrix1;
    if (matrix1.Col == matrix2.Col && matrix1.Row == matrix2.Row)
    {
#pragma omp parallel for
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

// friend function dotpow matrix dot power elementwise
matrix dotpow(const matrix &matrix, const double val)
{
    for (int i = 0; i < matrix.Numel; i++)
        matrix.Elem[i] = pow(matrix.Elem[i], val);
    return matrix;
}

// #######################
// friend operator overriding
// #######################

// friend operator + matrix add
matrix operator+(const matrix &matrix1, const matrix &matrix2)
{
    matrix temp(matrix1.Row, matrix1.Col);
    temp = matrix1;
    if (matrix1.Col == matrix2.Col && matrix1.Row == matrix2.Row)
    {
#pragma omp parallel for
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
    matrix Temp(matrix1.Row, matrix1.Col);
    Temp = matrix1;
#pragma omp parallel for
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
    matrix temp(matrix1.Row, matrix1.Col);
    temp = matrix1;
    if (matrix1.Col == matrix2.Col && matrix1.Row == matrix2.Row)
    {
#pragma omp parallel for
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
    matrix Temp(matrix1.Row, matrix1.Col);
    Temp = matrix1;
#pragma omp parallel for
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
    matrix Temp(matrix1.Row, matrix2.Col);
    double temp = 0.0;
    int j, k;
    if (matrix1.Col == matrix2.Row)
    {
#pragma omp parallel for private(temp, j, k)
        for (int i = 0; i < Temp.Row; i++)
        {
            for (j = 0; j < Temp.Col; j++)
            {
                temp = 0.0;
                for (k = 0; k < matrix1.Col; k++)
                {
                    // temp += matrix1.GetElem(i+1, k+1)*matrix2.Elem[k*matrix2.Col+j];
                    temp += matrix1.Elem[i * matrix1.Col + k] * matrix2.Elem[k * matrix2.Col + j];
                }
                Temp.Elem[i * Temp.Col + j] = temp;
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
    matrix Temp(matrix1.Row, matrix1.Col);
    Temp = matrix1;
#pragma omp parallel for
    for (int i = 0; i < matrix1.Numel; i++)
        Temp.Elem[i] = Temp.Elem[i] * val;
    return Temp;
}

// friend operator * value multiply (overloading 2)
matrix operator*(const double val, const matrix &matrix1)
{
    return matrix1 * val;
}

// friend operator / matrix dot divide elementwise
matrix operator/(const matrix &matrix1, const matrix &matrix2)
{
    matrix temp(matrix1.Row, matrix1.Col);
    temp = matrix1;
    if (matrix1.Col == matrix2.Col && matrix1.Row == matrix2.Row)
    {
        for (int i = 0; i < matrix1.Numel; i++)
            temp.Elem[i] /= matrix2.Elem[i];
        return temp;
    }
    else
    {
        printf("Error! Dimension is not match!");
        exit(1);
    }
}

// friend operator / value divide
matrix operator/(const matrix &matrix1, const double val)
{
    for (int i = 0; i < matrix1.Numel; i++)
        matrix1.Elem[i] = matrix1.Elem[i] / val;
    return matrix1;
}

#endif // __MATRIX_H__