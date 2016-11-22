#include <iostream>
#include <vector>
#include <cmath>

#include "fcn.hpp"
#include "composite_int.cpp"
#include "adaptive_int.cpp"

class fcn : public Fcn
{
public:
	double c, d;

	//Function evaluation.
	double operator()(double x) 
	{   
		return (exp(c*x) + sin(d*x));
	}

	//Antiderivate evaluation.
	double antiderivative(double x) 
	{
		return (exp(c*x)/c - cos(d*x)/d);
	}
};

//Tests Gauss method with 4 nodes on the function given above.
int main()
{
	//Limits of integration.
	double a = -3.0;
	double b = 5.0;

	//Integrand.
	fcn f;
	f.c = 0.5;
	f.d = 25.0;

	//True integral value.
	double If = f.antiderivative(b) - f.antiderivative(a);

	std::vector<double> rtol = {pow(10,-2), pow(10,-4), pow(10,-6), pow(10,-8),
		pow(10,-10), pow(10,-12)};
	
	std::vector<double> atol(rtol.size());
	for (int i = 0; i < atol.size(); i++)
		atol[i] = rtol[i] / 1000.0;
	
	double Rf, err, tol;
	int n, Ntot;

	std::cout << "\n  --------------------------------------\n\n";

	for (int i = 0; i < rtol.size(); i++)
	{
		adaptive_int(f, a, b, rtol[i], atol[i], Rf, n, Ntot);

		err = fabs(If - Rf);
		tol = rtol[i] * fabs(If) + atol[i];

		std::cout << "n: " << n << std::endl;
		std::cout << "Number of n's: " << Ntot << std::endl;
		std::cout << "Error: " << err << std::endl;
		std::cout << "Tolerance: " << tol << std::endl << std::endl;
	}

	std::cout << "  --------------------------------------\n\n";
}