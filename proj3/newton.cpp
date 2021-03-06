//*****************************************************************************
//	Author: Ben Maxey
//
//	The three functions in newton.cpp interpolate a given function using
//	a polynomial in Newton form.
//*****************************************************************************

#include "newton.hpp"

//Evaluates Newton basis function with nodes xnodes, degree n, at x
double Newton_basis(Matrix& xnodes, int n, double x)
{
	double product = 1.0;

	for (int i = 0; i <= n; i++)
		product *= x - xnodes(i);

	return product;
}

//Evaluates Newton polynomial with coefficients a, nodes xnodes, at x
double Newton_nestedform(Matrix& a, Matrix& xnodes, double x)
{
	int n = a.Cols() - 1;
	double output = a(n);

	for (int i = n - 1; i >= 0; i--)
		output = a(i) + output * (x - xnodes(i));

	return output;
}

//Generates coefficients iteratively using previously calculated coefficients
Matrix Newton_coefficients(Matrix& xnodes, Matrix& ynodes)
{
	Matrix a(1,xnodes.Cols());
	a(0) = ynodes(0);

	for (int i = 1; i < xnodes.Cols(); i++)
	{
		double pnx = Newton_nestedform(a, xnodes, xnodes(i));
		double basis = Newton_basis(xnodes, i - 1, xnodes(i));
		a(i) = (ynodes(i) - pnx) / basis;
	}

	return a;
}