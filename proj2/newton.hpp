//*****************************************************************************
//	Author:	Ben Maxey
//
//	This method solves a nonlinear equation using Newton's method.  The
//	method takes as arguments a function and its derivative, an initial guess,
//  a max number of iterations, an error tolerance, and a flag that determines
//	whether output should be printed.
//*****************************************************************************

#ifndef NEWTON_H
#define NEWTON_H

#include <cmath>
#include <iostream>

#include "fcn.hpp"

double newton(Fcn& f, Fcn& df, double x, int maxit, 
			  double tol, bool show_iterates);

#endif