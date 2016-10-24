#include "Lagrange2D.h"
#include "Lagrange.cpp"

double Lagrange2D(Matrix& x, Matrix& y, Matrix& f, double a, double b)
{
	double output = 0.0;

	int m = x.Cols() - 1;
	int n = y.Cols() - 1;

	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double xbasis = Lagrange_basis(x, m, a);
			double ybasis = Lagrange_basis(y, n, b);
			output += xbasis * ybasis * f(a, b);
		}
	}

	return output;
}