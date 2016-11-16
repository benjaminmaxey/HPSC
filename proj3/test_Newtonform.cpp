//*****************************************************************************
//	Author: Ben Maxey
//
//	test_Newtonform tests the functions in newton.cpp by using them to
//	interpolate a fourth degree polynomial.
//*****************************************************************************

#include <iostream>
#include <cmath>

#include "fcn.hpp"
#include "newton.hpp"

//Define test polynomial
class testFcn : public Fcn
{
public:
	double operator() (double x)
	{
		return 3.1 * pow(x, 4) + 2.3 * pow(x, 3) - 6.6 * pow(x, 2) +
			8.7 * x + 7.9;
	}
};

int main(int argc, char** argv)
{
	testFcn f;

	//Generate nodes
	double xvalues[5] = {-2, -1, 0, 1, 2};
	Matrix xnodes(1, 5, xvalues);
	Matrix ynodes(1, 5);

	//Get function value at each node
	for (int i = 0; i < 5; i++)
		ynodes(i) = f(xnodes(i));

	//Generate coefficients of Newton polynomial
	Matrix a = Newton_coefficients(xnodes, ynodes);

	Matrix x = Linspace(-3, 3, 1, 201);
	Matrix fx(1, 201);
	Matrix pnx(1, 201);
	Matrix err(1, 201);

	//Evaluate all polynomials and determine error
	for (int i = 0; i < 201; i++)
	{
		fx(i) = f(x(i));
		pnx(i) = Newton_nestedform(a, xnodes, x(i));
		err(i) = fx(i) - pnx(i);
	}

	//Write results to disk
	x.Write("x.txt");
	fx.Write("fx.txt");
	pnx.Write("pnx.txt");
	err.Write("err.txt");
}