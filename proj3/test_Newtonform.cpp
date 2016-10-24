#include <iostream>
#include <cmath>

#include "fcn.hpp"
#include "newton.hpp"

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

	double xvalues[5] = {-2, -1, 0, 1, 2};
	Matrix xnodes(1, 5, xvalues);
	Matrix ynodes(1, 5);

	for (int i = 0; i < 5; i++)
		ynodes(i) = f(xnodes(i));

	Matrix a = Newton_coefficients(xnodes, ynodes);

	Matrix x = Linspace(-3, 3, 1, 201);
	Matrix fx(1, 201);
	Matrix pnx(1, 201);
	Matrix err(1, 201);

	for (int i = 0; i < 201; i++)
	{
		fx(i) = f(x(i));
		pnx(i) = Newton_nestedform(a, xnodes, x(i));
		err(i) = fx(i) - pnx(i);
	}

	x.Write("x.txt");
	fx.Write("fx.txt");
	pnx.Write("pnx.txt");
	err.Write("err.txt");
}