//*****************************************************************************
//	Author: Ben Maxey
//
//	The three functions in newton.cpp interpolate a given function using
//	a polynomial in Newton form.
//*****************************************************************************

#ifndef NEWTON_H
#define NEWTON_H

#include "matrix.hpp"

double Newton_basis(Matrix& xnodes, int n, double x);
double Newton_nestedform(Matrix& a, Matrix& xnodes, double x);
Matrix Newton_coefficients(Matrix& xnodes, Matrix& ynodes);

#endif