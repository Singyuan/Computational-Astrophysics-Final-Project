#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <iostream>
#include <math.h>
class matrix
{
private:
    int Col;           //numeber of column of matrix
    int Row;           //numeber of row of matrix
    int Numel;         //numeber of element of matrix
    double *Elem;      //element of matrix
public:
    matrix(void);                                    //contruct matrix by default
    matrix(int Rows, int Cols, double Val = 0);      //contruct matrix by particular numel and default value 0
    ~matrix(void);

    void display();                                  //display matrix
    double GetElem(int R, int C) const;             //view the value of (R,C)
    void SetElem(int R, int C, double val);         //set the value of (R,C)
    matrix T();                                // transpose

    double &operator()(int R, int C);
    matrix &operator=(const matrix &matrix);   // given

public:
    friend matrix dotprod(const matrix &matrix1, const matrix &matrix2);     // matrix dot product elementwise
    friend matrix dotpow(const matrix &matrix, const double val);            // matrix dot power elementwise

    friend matrix operator+(const matrix &matrix1, const matrix &matrix2);   // matrix add
    friend matrix operator+(const matrix &matrix1, const double val);        // value add (overloading 1)
    friend matrix operator+(const double val, const matrix &matrix1);        // value add (overloading 2)

    friend matrix operator-(const matrix &matrix1, const matrix &matrix2);   // matrix minus
    friend matrix operator-(const matrix &matrix1, const double val);        // value minus (overloading 1)
    friend matrix operator-(const double val, const matrix &matrix1);        // value minus (overloading 1)

    friend matrix operator*(const matrix &matrix1, const matrix &matrix2);   // matrix multiply LA
    friend matrix operator*(const matrix &matrix1, const double val);        // value multiply (overloading 1)
    friend matrix operator*(const double val, const matrix &matrix1);        // value multiply (overloading 2)
    friend matrix operator*(const matrix &matrix1, const matrix &matrix2);

    friend matrix operator/(const matrix &matrix1, const matrix &matrix2);   // matrix dot divide elementwise
    friend matrix operator/(const matrix &matrix1, const double val);        // value divide
    // friend std::ostream &operator<<(std::ostream &os, const matrix &matrix);
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
matrix  matrix::T()
{
	matrix Temp(Col, Row);
	for(int i = 0; i < Temp.Row; i++)
		for(int j = 0; j < Temp.Col; j++)
			Temp.Elem[i * Temp.Col +j] = Elem[j * Temp.Row +i];
	return Temp;
}

// display
void matrix::display(void)
{
    for (int i = 1; i <= Row; i++)
    {
        for (int j = 1; j <= Col; j++)
            printf("%f ", (*this)(i,j));
        printf("\n");
    }
}

// view the value of (R,C)
double matrix::GetElem(int R, int C) const
{
    return Elem[(R-1)*Col+(C-1)];
}
// set the value of (R,C)
void matrix::SetElem(int R, int C, double val)
{
    Elem[(R - 1) * Col + (C - 1)] = val;
}



// #######################
// operator overriding
// #######################

// operator =
matrix &matrix::operator=(const matrix &M)
{
    if (this == &M)
    {
        return *this;
    }
    Col = M.Col;
    Row = M.Row;
    Numel = M.Numel;
    Elem = new double[Numel];
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
        for (int i = 0; i < matrix1.Numel; i++)
            temp.Elem[i] *= matrix2.Elem[i];
        return temp;
    }
    else
    {
        printf("Error! Dimension is not match!");
        return temp;
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
        for (int i = 0; i < matrix1.Numel; i++)
            temp.Elem[i] += matrix2.Elem[i];
        return temp;
    }
    else
    {
        printf("Error! Dimension is not match!");
        return temp;
    }
}

//friend operator + value add (overloading 1)
matrix operator+(const matrix &matrix, const double val)
{
    for (int i = 0; i < matrix.Numel; i++)
        matrix.Elem[i] += val;
    return matrix;
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
        for (int i = 0; i < matrix1.Numel; i++)
            temp.Elem[i] -= matrix2.Elem[i];
        return temp;
    }
    else
    {
        printf("Error! Dimension is not match!");
        return temp;
    }
}

//friend operator - value minus (overloading 1)
matrix operator-(const matrix &matrix, const double val)
{
    for (int i = 0; i < matrix.Numel; i++)
        matrix.Elem[i] -= val;
    return matrix;
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
    if (matrix1.Col == matrix2.Row)
    {
        for (int i = 0; i < Temp.Row; i++)
            for (int j = 0; j < Temp.Col; j++)
            {
                double temp = 0.0;
                for (int k = 0; k < matrix1.Col; k++)
                {
                    temp += matrix1.GetElem(i+1, k+1)*matrix2.Elem[k*matrix2.Col+j];
                }
                Temp.SetElem(i + 1, j + 1, temp);
            }
    return Temp;
    }
    else
    {
        printf("Error! Dimension is not match!");
        return matrix1;
    }
}

//friend operator * value multiply (overloading 1)
matrix operator*(const matrix &matrix1, const double val)
{
    for (int i = 0; i < matrix1.Numel; i++)
        matrix1.Elem[i] = matrix1.Elem[i] * val;
    return matrix1;
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
        return temp;
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