//*****************************************************************************
//	Author:	Ben Maxey
//
//	This method solves a nonlinear equation using Newton's method.  The
//	method takes as arguments a function and its derivative, an initial guess,
//  a max number of iterations, an error tolerance, and a flag that determines
//	whether output should be printed.
//*****************************************************************************

#include "newton.hpp"

double newton(Fcn& f, Fcn& df, double x, int maxit, 
			  double tol, bool show_iterates)
{	
	if (show_iterates)
	{
		std::cout << "\n----------------------------------" << std::endl
				  << "----------------------------------" << std::endl
				  << "Initial guess: " << x << std::endl
				  << "Maximum iterations: " << maxit << std::endl
				  << "Tolerance: " << tol << std::endl
				  << "----------------------------------" << std::endl
				  << "----------------------------------\n" << std::endl;
	}

	double fx = f(x);
	
	for (int i = 0; i < maxit; i++)
	{
		double fp = df(x);
		double d = fx/fp;
		x -= d;
		fx = f(x);

		if (d < 0)
			d = 0 - d;

		if (show_iterates)
		{
			std::cout << "Iteration " << i + 1 << ": " << std::endl
					  << "\tCurrent guess: x = " << x << std::endl
					  << "\tSolution update: |h| = " << d << std::endl
					  << "\tCurrent residual: f(x) = " << fx << std::endl
					  << std::endl;
		}

		if (d < tol)
		{
			if (show_iterates)
				std::cout << "\tConverged at x = " << x << std::endl;
			
			return x;
		}
	}
}