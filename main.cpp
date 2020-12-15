#include <iostream>
#include <iomanip>
#include <random>
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

	size_t iterations = 5;

	int dim_u = 10;
	int dim_y = 1;
	int dim_x = 2;

	MatrixXd A(dim_x,dim_x);
	MatrixXd B(dim_x,dim_y);

	MatrixXd C(dim_y,dim_x);
	MatrixXd D(dim_y,dim_u);

VectorXd next_u(MatrixXd &x, VectorXd w)
{
	MatrixXd u(dim_u,1);

	for (int i = 0; i < dim_u; i++)
	{
		x.col(i) = A*x.col(i) + B*w(i);
		u.row(i) = C*x.col(i);
	}

	return u;
}

int main()
{
	A << 0.5, 0.5,
		  0.5, 0.5;

	B << 0.5,
		  0.5;

	C << 0.5, 0.5;

	VectorXd ug(dim_u);
	VectorXd uh(dim_u);

	MatrixXd xg(dim_x,dim_u);
	MatrixXd xh(dim_x,dim_u);

	VectorXd y(dim_u);

	// White Noise vector;
	default_random_engine generator;
	normal_distribution<double>* distribution;
		distribution = new normal_distribution<double>[2*dim_u];

	VectorXd wg(dim_u);
	VectorXd wh(dim_u);
		for (int n = 0; n < dim_u; ++n)
		{
			wg(n) = distribution[n](generator);
			wh(n) = distribution[n + dim_u](generator);
		}

// Correlating noise
	for (int k = 0; k < iterations; ++k)
	{
		ug = next_u(xg, wg);
		uh = next_u(xh, wh);

		/*
		cout << ug << endl << endl;
		cout << uh << endl << endl;
		*/

		for (int n = 0; n < dim_u; ++n)
		{
			wg(n) = distribution[n](generator);
			wh(n) = distribution[n + dim_u](generator);
		}

	}
	// How to make structure to hold x's in it?





	return 0;
}
