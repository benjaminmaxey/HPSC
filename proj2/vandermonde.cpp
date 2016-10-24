//*****************************************************************************
//	Author: Ben Maxey
//
//	This program solves a linear system Ax = b, where A is an (n x n)
//	Vandermonde matrix.  The system is solved using Gaussian elimination with
//	partial pivoting.  Values of n used are 5, 9, 17, 33, and 65.
//*****************************************************************************

#include <cmath>

#include "matrix.hpp"

//*****************************************************************************
//	Function definitions.
//*****************************************************************************

//Generates a Vandermonde matrix of order (size).
Matrix generate(size_t size)
{
	Matrix v = Linspace(0, 1, 1, size);
	Matrix A(size, size);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
			A(i, j) = pow(v(i), j);
	}

	return A;
}

//Solves the system Ax = b for x.  Passes arguments by value so that the
//original copies remain unaltered and can be used later for error calculation.
Matrix solve(Matrix A, Matrix b)
{
	Matrix x(b.Rows(), 1);

	//Gaussian elimination with partial pivoting.
	for (int k = 0; k < A.Cols() - 1; k++)
	{
		//Find largest entry in current column.
		int m = k;
		for (int i = k + 1; i < A.Cols(); i++)
		{
			if (abs(A(i, k) > abs(A(m,k))))
				m = i;
		}

		//Swap row containing max entry with row k.
		for (int i = k; i < A.Cols(); i++)
		{
			double temp = A(m, i);
			A(m, i) = A(k, i);
			A(k, i) = temp;
		}

		double temp = b(m, 0);
		b(m, 0) = b(k, 0);
		b(k, 0) = temp;

		//Carry out elimination.
		for (int i = k + 1; i < A.Cols(); i++)
		{
			A(i, k) /= A(k, k);

			for (int j = k + 1; j < A.Cols(); j++)
				A(i, j) -= A(i, k) * A(k, j);

			b(i, 0) -= A(i, k) * b(k, 0);
		}
	}

	//Backward substitution.
	for (int i = A.Cols() - 1; i >= 0; i--)
	{
		double sum = 0;

		for (int k = i + 1; k < A.Cols(); k++)
			sum += A(i, k) * x(k, 0);

		x(i, 0) = (b(i, 0) - sum) / A(i, i);
	}

	return x;
}

//Generates all matrices and vectors, solves, then calculates error and
//residual vectors.
int main(int argc, char* argv[])
{
	Matrix A5 = generate(5);
	Matrix A9 = generate(9);
	Matrix A17 = generate(17);
	Matrix A33 = generate(33);
	Matrix A65 = generate(65);

	Matrix x5 = Random(5, 1);
	Matrix x9 = Random(9, 1);
	Matrix x17 = Random(17, 1);
	Matrix x33 = Random(33, 1);
	Matrix x65 = Random(65, 1);

	Matrix b5 = A5 * x5;
	Matrix b9 = A9 * x9;
	Matrix b17 = A17 * x17;
	Matrix b33 = A33 * x33;
	Matrix b65 = A65 * x65;

	Matrix sol5 = solve(A5, b5);
	Matrix sol9 = solve(A9, b9);
	Matrix sol17 = solve(A17, b17);
	Matrix sol33 = solve(A33, b33);
	Matrix sol65 = solve(A65, b65);

	Matrix err5 = x5 - sol5;
	Matrix err9 = x9 - sol9;
	Matrix err17 = x17 - sol17;
	Matrix err33 = x33 - sol33;
	Matrix err65 = x65 - sol65;

	Matrix res5 = A5 * err5;
	Matrix res9 = A9 * err9;
	Matrix res17 = A17 * err17;
	Matrix res33 = A33 * err33;
	Matrix res65 = A65 * err65;

	std::cout << "Error norm with n = 5: " << Norm(err5) << std::endl;
	std::cout << "Error norm with n = 9: " << Norm(err9) << std::endl;
	std::cout << "Error norm with n = 17: " << Norm(err17) << std::endl;
	std::cout << "Error norm with n = 33: " << Norm(err33) << std::endl;
	std::cout << "Error norm with n = 65: " << Norm(err65) << std::endl;

	std::cout << "Resid norm with n = 5: " << Norm(res5) << std::endl;
	std::cout << "Resid norm with n = 9: " << Norm(res9) << std::endl;
	std::cout << "Resid norm with n = 17: " << Norm(res17) << std::endl;
	std::cout << "Resid norm with n = 33: " << Norm(res33) << std::endl;
	std::cout << "Resid norm with n = 65: " << Norm(res65) << std::endl;

	Matrix A5i = Inverse(A5);
	Matrix A9i = Inverse(A9);
	Matrix A17i = Inverse(A17);
	//A33 and A65 are singular; cannot determine inverses

	double k5 = Norm(A5) * Norm(A5i);
	double k9 = Norm(A9) * Norm(A9i);
	double k17 = Norm(A17) * Norm(A17i);

	std::cout << "Condition number with n = 5: " << k5 << std::endl;
	std::cout << "Condition number with n = 9: " << k9 << std::endl;
	std::cout << "Condition number with n = 17: " << k17 << std::endl;
}