//*****************************************************************************
//	Author: Ben Maxey
//
//	This program tests the functionality of the newton() function.  The
//	function to be solved is x^2(x - 3)(x + 2).  The initial guesses that will
//	be used are x = -3, 1, 2, and the tolerances are 10^(-1), 10^(-5), and 
//	10^(-9).
//*****************************************************************************

#include <iostream>
#include <cmath>

#include "newton.hpp"
#include "fcn.hpp"

//*****************************************************************************
//	Subclass definitions.
//*****************************************************************************

//f(x) = x^2(x - 3)(x + 2) = x^4 - x^3 - 6x^2
class Func : public Fcn
{
public:
	double operator() (double x)
	{
		return pow(x, 2) * (x - 3) * (x + 2);
	}
};

//f'(x) = 4x^3 - 3x^2 - 12x
class Deriv : public Fcn
{
public:
	double operator() (double x)
	{
		return 4 * pow(x, 3) - 3 * pow(x, 2) - 12 * x;
	}
};

//*****************************************************************************
//	Function definition.
//*****************************************************************************

int main(int argc, char* argv[])
{
	Func fx;
	Deriv fp;	

	double z1 = newton(fx, fp, -3, 50, pow(10, -1), true);
	double z2 = newton(fx, fp, 1, 50, pow(10, -1), true);
	double z3 = newton(fx, fp, 2, 50, pow(10, -1), true);
	double z4 = newton(fx, fp, -3, 50, pow(10, -5), true);
	double z5 = newton(fx, fp, 1, 50, pow(10, -5), true);
	double z6 = newton(fx, fp, 2, 50, pow(10, -5), true);
	double z7 = newton(fx, fp, -3, 50, pow(10, -9), true);
	double z8 = newton(fx, fp, 1, 50, pow(10, -9), true);
	double z9 = newton(fx, fp, 2, 50, pow(10, -9), true);
}