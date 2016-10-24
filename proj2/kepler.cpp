//*****************************************************************************
//	Author: Ben Maxey
//
//	This program models the motion of an object in an elliptical orbit.  The
//	major radius of the ellipse is 2 while the minor radius is 1.25.  The angle
//	of the object around its orbit is determined for t = [0,10] using Kepler's
//	equation, then the position of the object is determined using the equation
//	for the radial position of an object in an elliptical orbit.
//*****************************************************************************

#include <iostream>
#include <cmath>

#include "matrix.hpp"
#include "newton.hpp"
#include "fcn.hpp"

//*****************************************************************************
//	Subclass definitions.
//*****************************************************************************

//f(x) = [(1 - b^2/a^2)^(1/2)]sin(x) - x - t
//b = 1.25, a = 2, t = time
class Func : public Fcn
{
private:
	double t;
	double e;
public:
	Func()
	{
		t = 0;
		e = sqrt(1 - (1.25 * 1.25)/(4));
	}

	Func(double time)
	{
		t = time;
		e = sqrt(1 - (1.25 * 1.25)/(4));
	}

	double operator() (double x)
	{
		return e * sin(x) - x - t;
	}
};

//f'(x) = [(1 - b^2/a^2)^(1/2)]cos(x) - 1
class Deriv : public Fcn
{
private:
	double e = sqrt(1 - (1.25 * 1.25)/(4));
public:
	double operator() (double x)
	{
		return e * cos(x) - 1;
	}
};

//r(x) = ab/{[(bcos(x))^2 + (asin(x))^2]^(1/2)}
//b = 1.25, a = 2
class Radial : public Fcn
{
public:
	double operator() (double x)
	{
		double denominator = sqrt(pow(1.25 * cos(x), 2) + pow(2 * sin(x), 2));
		return 2.5 / denominator;
	}
};

//*****************************************************************************
//	Function definitions.
//*****************************************************************************

//Uses newton() function to determine a zero of Kepler's equation for each of
//the given times.  The initial guess for each solve is the result of the
//previous solve.
Matrix solve(int initialGuess, Matrix& times)
{
	Matrix output(1, times.Cols());
	double initial = initialGuess;

	for (int i = 0; i < times.Cols(); i++)
	{
		double time = times(i);

		Func f(time);
		Deriv df;

		output(i) = initial = newton(f, df, initial, 6, pow(10, -5), false);
	}

	return output;
}

//Calculates the x value of the object at each angle.
Matrix computeX(Matrix& angles)
{
	Matrix output(1, angles.Cols());
	Radial r;

	for (int i = 0; i < angles.Cols(); i++)
		output(i) = r(angles(i)) * cos(angles(i));

	return output;
}

//Calculates the y value of the object at each angle.
Matrix computeY(Matrix& angles)
{
	Matrix output(1, angles.Cols());
	Radial r;

	for (int i = 0; i < angles.Cols(); i++)
		output(i) = r(angles(i)) * sin(angles(i));

	return output;
}

//Generates times from 0 to 10, then determines angles and positions.
int main(int argc, char** argv)
{
	Matrix times = Linspace(0, 10, 1, pow(10,4) + 1);
	times.Write("t.txt");

	Matrix w = solve(0, times);
	w.Write("w.txt");

	Matrix x = computeX(w);
	x.Write("x.txt");

	Matrix y = computeY(w);
	y.Write("y.txt");

	std::cout << "Generated files: t.txt w.txt x.txt y.txt" << std::endl;
}