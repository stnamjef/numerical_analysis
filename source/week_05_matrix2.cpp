#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

double Loss(double* y, double* y_hat, int N)
{
	double sse = 0.0;
	for (int i = 0; i < N; i++) {
		sse += pow((y[i] - y_hat[i]), 2);
	}
	return sse;
}

int main() {
	// read data
	ifstream data("data1.txt");
	int N;
	data >> N;

	double* x = new double[N];
	double* y = new double[N];

	for (int i = 0; i < N; i++) {
		data >> x[i] >> y[i];
	}
	data.close();

	// LSM1: using two equations
	MatrixXd A1(2, 2);
	VectorXd b1(2), c1(2);

	double Nx, Nxx, Ny, Nxy;
	Nx = Nxx = Ny = Nxy = 0;
	for (int i = 0; i < N; i++) {
		Nx += x[i];
		Ny += y[i];
		Nxx += x[i] * x[i];
		Nxy += x[i] * y[i];
	}

	A1 << N, Nx, Nx, Nxx;
	b1 << Ny, Nxy;

	c1 = A1.inverse() * b1;

	// LSM2: A * c = b -> c = (At*A)^(-1) * At * b
	MatrixXd A2(N, 2);
	VectorXd b2(N), c2(2);

	for (int i = 0; i < N; i++) {
		A2(i, 0) = 1;
		A2(i, 1) = x[i];
		b2[i] = y[i];
	}

	MatrixXd At = A2.transpose();

	c2 = (At * A2).inverse() * At * b2;

	// LSM3: gradient descent
	double* yhat = new double[N];

	double lr = 0.0000001, eta = 0.001;
	double a0 = 0.0, a1 = 0.0;
	double da0 = 0.0, da1 = 0.0;

	int iter = 0; double loss = 1;
	while (loss > eta && iter < 1000) {
		cout << loss << endl;
		// calc yhat
		for (int i = 0; i < N; i++) {
			yhat[i] = a0 + a1 * x[i];
		}
		// calc loss
		loss = Loss(y, yhat, N);
		// calc derivatives of Loss() w.r.t a0, a1
		for (int i = 0; i < N; i++) {
			da0 -= 2 * (y[i] - yhat[i]);
			da1 -= 2 * x[i] * (y[i] - yhat[i]);
		}
		// update weights
		a0 -= lr * da0;
		a1 -= lr * da1;
		// reset derivatives
		da0 = da1 = 0.0;
		// update iter
		iter++;
	}

	// Print results
	cout << "LSM1: " << c1.transpose() << endl;
	cout << "LSM2: " << c2.transpose() << endl;
	cout << "LSM3: " << a0 << ' ' << a1 << endl;

	delete[] x;
	delete[] y;
	delete[] yhat;
}