#include <iostream>
#include <vector>
#include <cmath>

#include "fcn.hpp"
#include "composite_int.cpp"

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
	printf("\n True Integral = %22.16e\n", If);

	//Test Gauss method with 4 nodes.
	std::cout << "\n Gauss-4 Approximation:\n";
	std::cout << "     n             R(f)            RelErr    Conv Rate\n";
	std::cout << "  ---------------------------------------------------\n";
 
	std::vector<int> n = {40, 80, 120, 160, 200, 240, 280, 320, 360};
	std::vector<double> errors(n.size());
	std::vector<double> hvals(n.size());

	//Iterate over n values and print approximations, errors, and
	//convergence rates.
	double Rf;
	for (int i = 0; i < n.size(); i++)
	{
		printf("   %6i", n[i]);

		Rf = composite_int(f, a, b, n[i]);
		errors[i] = fabs(If - Rf)/fabs(If);
		hvals[i] = (b - a)/n[i];
		

		if (i == 0) 
			printf("  %22.16e  %7.1e     ----\n", Rf, errors[i]);
		else
		{
			double convergence = (log(errors[i - 1]) - log(errors[i])) /
				(log(hvals[i - 1]) - log(hvals[i]));
			printf("  %22.16e  %7.1e   %f\n", Rf, errors[i], convergence);
		} 
	}

	std::cout << "  ---------------------------------------------------" <<
		std::endl;
}