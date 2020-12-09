#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
using namespace std;

void set_zero(double* ary, int n) { for (int i = 0; i < n; i++) ary[i] = 0.; }

double binomial(double n, double x, double p)
{
	double prob = std::pow(p, x) * std::pow((1. - p), (n - x));
	double upper = (x > n - x) ? x : n - x;
	double i, j;
	for (i = upper + 1, j = 1.; i <= n; i++, j++) {
		prob *= i / j;
	}
	return prob;
}

double poisson(double alpha, double x)
{
	double p = 1.0;
	for (double i = 1.; i <= x; i++) {
		p *= alpha / i;
	}
	return p * std::exp(-alpha);
}

// Binomial distribution: N packets problem
int main()
{
	int n = 10, N = 100000;
	double p = 0.1;

	// Analytical solution
	double* Pk1 = new double[n + 1];
	set_zero(Pk1, n + 1);

	double Ek1 = n * p, Vk1 = n * p * (1 - p);
	for (int k = 0; k <= n; k++) {
		Pk1[k] = binomial(n, k, p);
	}

	// Numerical solution
	double* Pk2 = new double[n + 1];
	set_zero(Pk2, n + 1);
	
	// 여기서부터 다시
	double Ek2 = 0.0, Vk2 = 0.0;
	for (int i = 0; i < N; i++) {
		int k = 0;
		for (int j = 0; j < n; j++) {
			if (rand() / (double)RAND_MAX < p) {
				k++;
			}
		}
		Pk2[k] += 1;
		Ek2 += k;
		Vk2 += k * k;
	}

	for (int k = 0; k <= n; k++) {
		Pk2[k] /= (double)N;
	}

	Ek2 = Ek2 / N;
	Vk2 = Vk2 / N - Ek2 * Ek2;

	// Poisson approximation
	double* Pk3 = new double[n + 1];
	set_zero(Pk3, n + 1);

	double Ek3 = 0.0, Vk3 = 0.0;
	for (int k = 0; k <= n; k++) {
		Pk3[k] = poisson(n * p, k);
		Ek3 += k * Pk3[k];
		Vk3 += k * k * Pk3[k];
	}
	Vk3 -= Ek3 * Ek3;

	cout << std::fixed << setprecision(4);
	cout << "Probabilities:" << endl;
	for (int k = 0; k <= n; k++) {
		cout << Pk1[k] << ", " << Pk2[k] << ", " << Pk3[k] << endl;
	}
	cout << "Mean: " << Ek1 << ", " << Ek2 << ", " << Ek3 << endl;
	cout << "Var : " << Vk1 << ", " << Vk2 << ", " << Vk3 << endl;

	delete[] Pk1;
	delete[] Pk2;
	delete[] Pk3;

	return 0;
}


//// sum of 5 dices
//
//double gaussian(double mu, double sigma, double x)
//{
//	double inv_sqrt_2pi = 0.3989422804014327;
//	double a = (x - mu) / sigma;
//	return inv_sqrt_2pi / sigma * std::exp(-0.5f * a * a);
//}
//
//int main()
//{
//	int n = 10, N = 100000;
//	int lower = n, upper = 6 * n;
//	int size = upper - lower + 1;
//
//	double mean = 0., var = 0.;
//	double* P1 = new double[size];
//	set_zero(P1, size);
//
//	// random simulation
//	int* dice = new int[n];
//	for (int i = 0; i < N; i++) {
//		for (int j = 0; j < n; j++) {
//			dice[j] = rand() % 6 + 1;
//		}
//		int sum = 0;
//		for (int j = 0; j < n; j++) {
//			sum += dice[j];
//		}
//		P1[sum - n] += 1.;
//		mean += (double)sum;
//		var += (double)(sum * sum);
//	}
//
//	// calculate probabilities
//	for (int i = 0; i < size; i++) {
//		P1[i] /= (double)N;
//	}
//
//	// calculate mean & var
//	mean /= (double)N;
//	var = var / (double)N - mean * mean;
//
//	// Gaussian simulation
//	double* P2 = new double[size];
//
//	// sum of one dice -> mu = 3.5, var = 2.91
//	// sum of n dices -> mu = n * 3.5, var = sqrt(n) * 2.91
//	double mu = n * 3.5;
//	double sigma = std::sqrt(n * 2.91);
//	for (int i = 0; i < size; i++) {
//		// RV -> sum of n dices
//		// if n == 6, then 6, 7, ..., 60
//		P2[i] = gaussian(mu, sigma, i + n);
//	}
//
//	cout << std::fixed << std::setprecision(4);
//	for (int i = 0; i < size; i++) {
//		cout << i + n << ": " << P1[i] << ", " << P2[i] << endl;
//	}
//	cout << std::setprecision(2);
//	cout << "Mean: " << mean << ", " << mu << endl;
//	cout << "Var: " << var << ", " << sigma * sigma << endl;
//
//	delete[] P1;
//	delete[] P2;
//
//	return 0;
//}


//// X ~ B(10, 1 / 6) simulation
//
//double fac(double x)
//{
//	if (x == 0.0 || x == 1.0) {
//		return 1.0;
//	}
//	return x * fac(x - 1);
//}
//
//double binomial(int n, int x, double p)
//{
//	return pow(p, x) * pow((1. - p), (n - x)) * fac(n) / fac(x) / fac(n - x);
//}
//
//// binomial
//int main()
//{
//	// X ~ B(10, 1 / 6) simulation
//	int n = 10, N = 10000;
//	double p = 1 / 6.;
//	double mean = 0., var = 0.;
//	double* P = new double[n + 1];
//	set_zero(P, n + 1);
//	for (int i = 0; i < N; i++) {
//		int x = 0;
//		for (int j = 0; j < n; j++) {
//			if (rand() / (double)RAND_MAX < p) x++;
//		}
//		P[x] += 1.;
//		mean += x;
//		var += x * x;
//	}
//
//	// get probability
//	for (int x = 0; x < n + 1; x++) P[x] /= (double)N;
//	mean /= (double)N;
//	var /= (double)N;
//	var -= mean * mean;
//
//	cout << "Mean: " << mean << ", " << n * p << endl;
//	cout << "Var: " << var << ", " << n * p * (1. - p) << endl;
//	cout << "Probabilities:" << endl;
//	for (int x = 0; x < n + 1; x++) {
//		cout << "x=" << x << ": " << P[x] << ", " << binomial(n, x, p) << endl;
//	}
//
//	delete[] P;
//
//	return 0;
//}


//// pi
//int main()
//{
//	int N = 10000;
//	float x, y, hit = 0;
//
//	ofstream out("pi.csv");
//	for (int i = 0; i < N; i++) {
//		x = (rand() / (float)RAND_MAX - 0.5f) * 2;
//		y = (rand() / (float)RAND_MAX - 0.5f) * 2;
//		/*x = rand() / (float)RAND_MAX;
//		y = rand() / (float)RAND_MAX;*/
//		if (x * x + y * y < 1.f) {
//			hit += 1.f;
//			out << x << ',' << y << endl;
//		}
//	}
//	out.close();
//
//	cout << hit / (float)N * 4.f << endl;
//
//	return 0;
//}


//// tetris, cdf -> prob
//int main()
//{
//	float p[] = { 12.3f, 26.5f, 17.23f, 22.3f, 9.2f, 12.47f };
//	float cdf[6];
//	
//	cdf[0] = 12.3f;
//	for (int i = 1; i < 6; i++) {
//		cdf[i] = cdf[i - 1] + p[i];
//	}
//
//	float count[6];
//	for (int i = 0; i < 10000; i++) {
//		float rr = rand() / (float)RAND_MAX * 100.f;
//		if (rr < cdf[0]) count[0] += 1.f;
//		else if (rr < cdf[1]) count[1] += 1.f;
//		else if (rr < cdf[2]) count[2] += 1.f;
//		else if (rr < cdf[3]) count[3] += 1.f;
//		else if (rr < cdf[4]) count[4] += 1.f;
//		else count[5] += 1.f;
//	}
//
//	for (int i = 0; i < 6; i++) {
//		count[i] /= 100.f;
//		cout << ' ' << i << ' ' << count[i] << ' ' << p[i] << endl;
//	}
//
//	return 0;
//}