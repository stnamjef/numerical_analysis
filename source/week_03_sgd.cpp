#include <iostream>
#include <fstream>
#include <functional>
#include <random>
#include <cmath>
using namespace std;

double fx(double a, double x, double b) { return a * x + b; }

double square_dist(double a, double b, double x, double y)
{ 
	// dist btw (x_i, y_i) and ax - y + b = 0
	// dist = |ax_i - y_i| / sqrt(a^2 + 1)
	return (a * x - y) * (a * x - y) / (a * a + 1);
}

double SSE(double a, double b, double* x, double* y, int n)
{
	double sse = 0;
	for (int i = 0; i < n; i++) {
		sse += square_dist(a, b, x[i], y[i]);
	}
	return sse;
}

double EE(double x0, double y0, double x1, double y1)
{
	return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
}

void randomNormal(double* x, int n, double mu, double var)
{
	default_random_engine e(444);
	normal_distribution<double> dist(0, var);
	std::for_each(x, x + n, [&](double& elem) { elem += dist(e); });
}

double dLossda(double a, double da, double b, double* x, double* y, int n)
{
	return (SSE(a + da, b, x, y, n) - SSE(a, b, x, y, n)) / da;
}

double dLossdb(double a, double b, double db, double* x, double* y, int n)
{
	return (SSE(a, b + db, x, y, n) - SSE(a, b, x, y, n)) / db;
}

int main()
{
	int n = 100;
	double* X = new double[n];
	double* Y = new double[n];

	// generate evenly spaced x coordinates ranged from 0 to 30
	double x = 0;
	double delta = 30. / (double)n;
	for (int i = 0; i < n; i++) {
		X[i] = x;
		x += delta;
	}

	// calculate y coordinates y =  2 * x + 1
	for (int i = 0; i < n; i++) {
		Y[i] = fx(2, X[i], 1);
	}

	// add random noise to the y coordinates
	randomNormal(Y, n, 0, 10);

	// write the x, y coordinates as a csv file
	ofstream out("data_original.csv");
	
	for (int i = 0; i < n; i++) {
		out << X[i] << ',' << Y[i] << endl;
	}

	out.close();

	// optimize the loss function
	// target values are a = 2, b = 1

	double psi = 0.0001, eta = 0.00001, da = 0.001, db = 0.001;
	double a0 = 0, b0 = 0, a1 = 0.2, b1 = 0.2;

	int iter = 0;
	while (EE(a0, b0, a1, b1) > eta && iter < 1000) {
		a0 = a1;
		b0 = b1;
		a1 = a1 - psi * dLossda(a0, da, b0, X, Y, n);
		b1 = b1 - psi * dLossdb(a0, b0, db, X, Y, n);

		cout << "Iter: " << iter + 1 << ", Error: " << EE(a0, b0, a1, b1) << endl;

		iter++;
	}

	cout << "Final output a1: " << a1 << ", b1: " << b1 << endl;

	// estimate y_hat using the function derived above
	out.open("data_estimated.csv");

	for (int i = 0; i < n; i++) {
		out << X[i] << ',' << fx(a1, X[i], b1) << endl;
	}

	out.close();

	delete[] X;
	delete[] Y;

	return 0;
}