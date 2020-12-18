#include <iostream>
#include <iomanip>
#include <random>
#include "Eigen/Dense"

using namespace Eigen;
using namespace std;

	size_t iterations = 10000;

	int num_iput = 10;
	int dim_meas = 1;
	int dim_stat = 2;

	MatrixXd As(dim_stat,dim_stat);
	MatrixXd Bs(dim_stat,dim_meas);
	MatrixXd Cs(dim_meas,dim_stat);
	MatrixXd Ds(dim_meas,dim_meas);

	MatrixXd Ag(dim_stat,dim_stat);
	MatrixXd Bg(dim_stat,dim_meas);
	MatrixXd Cg(dim_meas,dim_stat);
	MatrixXd Dg(dim_meas,dim_meas);

	MatrixXd Ah(dim_stat,dim_stat);
	MatrixXd Bh(dim_stat,dim_meas);
	MatrixXd Ch(dim_meas,dim_stat);
	MatrixXd Dh(dim_meas,dim_meas);

VectorXd transfer(const MatrixXd A, const MatrixXd B, const MatrixXd C, const MatrixXd D, MatrixXd &state, VectorXd ivec)
{
	MatrixXd output(num_iput,1);

	for (int i = 0; i < num_iput; i++)
	{
		output.row(i) = C*state.col(i) + D*ivec(i);	// Computes current output vector component
		state.col(i) = A*state.col(i) + B*ivec(i);	// Increments state vector
	}

	return output;
}

int main()
{
	As << -0.5, 0.25,
		  0.25, -0.5;

	Bs << 0.5,
		  0.5;

	Cs << 0.5, 0.5;

	Ds << 0.0;

	Ag = Ah = As;
	Bg = Bh = Bs;
	Cg = Ch = Cs;
	Dg = Dh = Ds;

	VectorXd ug(num_iput);
	VectorXd uh(num_iput);

	MatrixXd xsg(dim_stat,num_iput);
	MatrixXd xsh(dim_stat,num_iput);

	MatrixXd xg(dim_stat,num_iput);
	MatrixXd xh(dim_stat,num_iput);

	size_t mic_num = num_iput;
	VectorXd yg(mic_num);
	VectorXd yh(mic_num);

	// White Noise vector;
	default_random_engine generator;
	normal_distribution<double>* distribution;
		distribution = new normal_distribution<double>[2*num_iput];

	VectorXd wg(num_iput);
	VectorXd wh(num_iput);
		for (int n = 0; n < num_iput; ++n)
		{
			wg(n) = distribution[n](generator);
			wh(n) = distribution[n + num_iput](generator);
		}

// Correlating noise
	for (int k = 0; k < iterations; ++k)
	{
		ug = transfer(As, Bs, Cs, Ds, xsg, wg);
		uh = transfer(As, Bs, Cs, Ds, xsh, wh);

		yg = transfer(Ag, Bg, Cg, Dg, xg, ug);
		yh = transfer(Ah, Bh, Ch, Dh, xh, uh);

		cout << yg + yh << endl << endl;

		for (int n = 0; n < num_iput; ++n)
		{
			wg(n) = distribution[n](generator);
			wh(n) = distribution[n + num_iput](generator);
		}

	}

	return 0;
}
