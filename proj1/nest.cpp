//*****************************************************************************
//	Author: Ben Maxey
//
//	Evaluates a polynomial at the given x value.  The Matrix a should contain
//	the coefficients of the polynomial in ascending order.
//*****************************************************************************

#include "nest.h"

double nest(Matrix& a, double x)
{
	int n = a.Cols() - 1;
	double output = a(n);

	for (int i = n - 1; i >= 0; i--)
		output = a(i) + (x * output);

	return output;
}