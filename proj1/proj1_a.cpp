//*****************************************************************************
//	Author: Ben Maxey
//	
//	Approximates the function e^x by generating three Taylor polynomials of
//	degree 4, 8, and 12.
//*****************************************************************************

#include <iostream>

#include "matrix.hpp"
#include "nest.h"

int factorial(int n);
Matrix exponentialTaylor(int degree);
Matrix evaluate(Matrix& x, Matrix& f);
Matrix error(Matrix& f, Matrix& p);

int main(int argc, char** argv)
{
	//Initialize matrix z.
	Matrix z = Linspace(-1.0, 1.0, 1, 201);
	z.Write("z.txt");

	//Initialize Taylor polynomial coefficient matrices.
	Matrix coeff4 = exponentialTaylor(4);
	Matrix coeff8 = exponentialTaylor(8);
	Matrix coeff12 = exponentialTaylor(12);

	//Evaluate coefficient matrices for all values in z.
	Matrix p4 = evaluate(z, coeff4);
	Matrix p8 = evaluate(z, coeff8);
	Matrix p12 = evaluate(z, coeff12);

	//Write p4, p8, and p12.
	p4.Write("p4.txt");
	p8.Write("p8.txt");
	p12.Write("p12.txt");

	//Evaluate exp(x) for all values in z, then write f.
	Matrix f(1, z.Cols());
	for (int i = 0; i < z.Cols(); i++)
		f(i) = exp(z(i));
	f.Write("f.txt");

	//Compute error for each vector and write the error vectors.
	Matrix err4 = error(f, p4);
	Matrix err8 = error(f, p8);
	Matrix err12 = error(f, p12);
	err4.Write("err4.txt");
	err8.Write("err8.txt");
	err12.Write("err12.txt");
}

//Returns the factorial of n.
int factorial(int n)
{
	if (n <= 0)
		return 1;

	int x = n, output = n;
	while (x > 1)
		output *= --x;

	return output;
}

//Returns a row vector with coefficients of the Taylor polynomial of the 
//specified degree for exp(x) at x = 0.
Matrix exponentialTaylor(int degree)
{
	Matrix output(1, degree);
	
	for (int i = 0; i < degree; i++)
		output(i) = 1.0/(factorial(i));

	return output;
}

//Returns a matrix whose entries are the value of f(x) for all entries in x,
//where f is a row vector containing the coefficients of a power series.
Matrix evaluate(Matrix& x, Matrix& f)
{
	Matrix output(1, x.Cols());

	for (int i = 0; i < x.Cols(); i++)
		output(i) = nest(f, x(i));

	return output;
}

//Returns a matrix whose entries are the absolute value of the difference
//between the values of f and p.
Matrix error(Matrix& f, Matrix& p)
{
	Matrix output(1, f.Cols());

	for (int i = 0; i < f.Cols(); i++)
	{
		if (f(i) - p(i) >= 0)
			output(i) = f(i) - p(i);
		else
			output(i) = p(i) - f(i);
	}
	
	return output;
}