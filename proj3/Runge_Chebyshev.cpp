#include <cmath>

#include "matrix.hpp"
#include "Lagrange2D.h"

class Runge
{
public:
	double operator() (double x, double y)
	{
		double denom = 1.0 + x * x + y * y;
		return 1.0/denom;
	}
};

class Chebyshev
{
public:
	double operator() (double L, double m, double i)
	{
		const double PI = 3.14159265;
		double argument = (2 * i + 1) * PI / (2 * m + 2);
		return L * cos(argument);
	}
};

int main(int argc, char** argv)
{
	Runge function;
	Chebyshev node;

	Matrix avals = Linspace(-4, 4, 1, 201);
	Matrix bvals = Linspace(-4, 4, 1, 101);

	Matrix x6(1,7);
	for (int i = 0; i < 7; i++)
		x6(i) = node(4, 6, i);

	Matrix y6(1,7);
	for (int i = 0; i < 7; i++)
		y6(i) = node(4, 6, i);

	Matrix f6(7,7);
	for (int i = 0; i < 7; i++)
	{
		for (int j = 0; j < 7; j++)
			f6(i,j) = function(x6(i), y6(j));
	}

	Matrix p6(201,101);
	for (int i = 0; i < 201; i++)
	{
		for (int j = 0; j < 101; j++)
			p6(i,j) = Lagrange2D(x6, y6, f6, avals(i), bvals(j));
	}

	Matrix x24(1,25);
	for (int i = 0; i < 25; i++)
		x24(i) = node(4, 24, i);

	Matrix y24(1,25);
	for (int i = 0; i < 25; i++)
		y24(i) = node(4, 24, i);

	Matrix f24(25,25);
	for (int i = 0; i < 25; i++)
	{
		for (int j = 0; j < 25; j++)
			f24(i,j) = function(x24(i), y24(j));
	}

	Matrix p24(201,101);
	for (int i = 0; i < 201; i++)
	{
		for (int j = 0; j < 101; j++)
			p24(i,j) = Lagrange2D(x24, y24, f24, avals(i), bvals(j));
	}

	p6.Write("p6.Cheb.txt");
	p24.Write("p24.Cheb.txt");
}