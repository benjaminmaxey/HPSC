#include "Lagrange2D.h"
#include "Lagrange.cpp"
#include <iostream>

double Lagrange2D(Matrix& x, Matrix& y, Matrix& f, double a, double b)
{
	double output = 0.0;

	int m = x.Size() - 1;
	int n = y.Size() - 1;

	for (int i = 0; i <= m; i++)
	{
		for (int j = 0; j <= n; j++)
		{
			double xbasis = Lagrange_basis(x, i, a);
			double ybasis = Lagrange_basis(y, j, b);
			output += xbasis * ybasis * f(i, j);
		}
	}

	return output;
}