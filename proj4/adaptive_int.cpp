int adaptive_int(Fcn& f, const double a, const double b, const double rtol,
	const double atol, double& R, int& n, int& Ntot)
{
	double Rn = composite_int(f, a, b, 1);
	double Rn_plus = composite_int(f, a, b, 2);

	Ntot = 2;
	int n_plus = 2;
	int n_minus = 2;
	
	double err, tol, step;
		
	while (Ntot < 100)
	{
		err = fabs(Rn_plus - Rn);
		tol = rtol * fabs(Rn_plus) + atol;

		if (err < tol)
		{
			R = Rn_plus;
			n = n_plus;
			return 0;
		}

		step = 0.5 * n_minus * pow((err/tol), 0.125);
		n_minus = n_plus;
		n_plus += step;

		Rn = Rn_plus;
		Rn_plus = composite_int(f, a, b, n_plus);
		Ntot++;
	}

	return 1;
}