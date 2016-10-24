//*****************************************************************************
//	Author: Ben Maxey
//
//	Estimates the derivative of x^3 at x=3 using a forward finite difference
//	approximation.  Also calculates the relative error and upper bound on the
//	relative error of these estimates.
//*****************************************************************************

#include <iostream>
#include <cmath>

#include "matrix.hpp"

Matrix estimate(Matrix& increments);
Matrix getRelError(Matrix& estimates);
Matrix getMaxError(Matrix& increments);

int main()
{
	Matrix n = Linspace(1, 52, 1, 52);
	n.Write("n.txt");

	//Generate increments
	Matrix h(1, 52);
	for (int i = 0; i < 52; i++)
		h(i) = exp2(-n(i));
	h.Write("h.txt");

	//Calculate forward difference approximations.
	Matrix estimates = estimate(h);
	estimates.Write("estimates.txt");

	//Calculate relative error.
	Matrix r = getRelError(estimates);
	r.Write("r.txt");

	//Calculate upper bound on relative error.
	Matrix R = getMaxError(h);
	R.Write("R.txt");
}

//Calculates the forward finite difference approximation for each given
//increments.
Matrix estimate(Matrix& increments)
{
	Matrix output(1, 52);

	for (int i = 0; i < 52; i++)
	{
		double z = 3 + increments(i);
		double fah = pow(z,-3);
		output(i) = (fah - (1.0/27.0)) / increments(i);
	}

	return output;
}

//Calculates the relative error for each estimate using the true value of
//d/dx(x^3) at x=3.
Matrix getRelError(Matrix& estimates)
{
	Matrix output(1, 52);
	double correct = -(1.0/27);

	for (int i = 0; i < 52; i++)
	{
		if (estimates(i) < correct)
			output(i) = correct - estimates(i);
		else
			output(i) = estimates(i) - correct;
	}

	return output;
}

//Calculates the upper bound on the relative error for each increment using the
//formula derived in the project handout.
Matrix getMaxError(Matrix& increments)
{
	Matrix output(1, 52);
	double c1 = (12.0 * pow(3.0,-5.0)) / (2.0 * (1.0/27.0));
	double c2 = (1.0/27.0) * pow(2,-52) / (1.0/27.0);

	for (int i = 0; i < 52; i++)
		output(i) = (c1 * increments(i)) + (c2 * 1.0/increments(i));

	return output;
}