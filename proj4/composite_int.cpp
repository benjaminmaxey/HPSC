//*****************************************************************************
//	Author: Ben Maxey
//
//	The composite_int function uses a composite Gaussian quadrature rule with
//	four nodes per subinterval to evaluate a function on a given interval.  The
//	function takes the following inputs: Fcn& f (the function to be
//	integrated), const double a, const double b (the lower and upper bounds of
//	integration), and const int n (the number of subintervals).
//*****************************************************************************

#include <iostream>
#include <cmath>

#include "fcn.hpp"

double composite_int(Fcn& f, const double a, const double b, const int n)
{
	//Width of intervals.
	double h = (b - a)/n;

	//Nodes and weights.
	double x0 = -sqrt((1/7) * (3 - 4 * sqrt(0.3)));
	double x1 = -sqrt((1/7) * (3 + 4 * sqrt(0.3)));
	double x2 = sqrt((1/7) * (3 - 4 * sqrt(0.3)));
	double x3 = sqrt((1/7) * (3 + 4 * sqrt(0.3)));

	double w0 = 0.5 + (1/12) * sqrt(10/3);
	double w1 = 0.5 - (1/12) * sqrt(10/3);

	double sum = 0.0;
	double xmid, n0, n1, n2, n3;

	//Apply Gaussian quadrature rule to each subinterval.
	for (int i = 0; i < n; i++)
	{
		//Determine evaluation points.
		xmid = a + (i + 0.5) * h
		n0 = xmid + 0.5 * h * x0;
		n1 = xmid + 0.5 * h * x1;
		n2 = xmid + 0.5 * h * x2;
		n3 = xmid + 0.5 * h * x3;

		sum += w0 * f(n0) + w1 * f(n1) + w0 * f(n2) + w1 * f(n3);
	}

	return sum * 0.5 * h;
}