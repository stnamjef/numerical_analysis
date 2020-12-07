#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
using namespace std;
using namespace Eigen;


double INF = 1000000000.;

void argsort(int* idx, int n, const EigenSolver<MatrixXd>& es)
{
	VectorXd data = es.eigenvalues().real();
	for (int i = 0; i < n; i++) {
		double max_value = -INF;
		int max_idx = -1;
		for (int j = 0; j < n; j++) {
			if (data[j] > max_value && data[j] != INF) {
				max_value = data[j];
				max_idx = j;
			}
		}
		data[max_idx] = INF;
		idx[i] = max_idx;
	}
}

void to_unsigned_char_range(MatrixXd& M)
{
	for (int i = 0; i < M.rows(); i++) {
		for (int j = 0; j < M.cols(); j++) {
			double value = M(i, j);
			if (value < 1e-10 || value > 255.) {
				M(i, j) = 0.;
			}
		}
	}
}

void save_matrix(string path, const MatrixXd& M)
{
	ofstream out(path);
	for (int i = 0; i < M.rows(); i++) {
		for (int j = 0; j < M.cols(); j++) {
			out << M(i, j) << ',';
		}
		out << endl;
	}
	out.close();
}

int main()
{
	MatrixXd A(7, 13);
	A << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
		0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0,
		0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0,
		0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
		0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
		0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0,
		0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;  // 7x13
	A *= 200;

	int m = A.rows();
	int s = A.cols();

	MatrixXd ATA = A.transpose() * A;

	EigenSolver<MatrixXd> es;
	es.compute(ATA);

	int n = es.eigenvalues().real().size();
	int* idx = new int[n];
	argsort(idx, n, es);

	MatrixXd V(s, s);
	for (int i = 0; i < s; i++) {
		V.col(i) = es.eigenvectors().col(idx[i]).real();
	}

	MatrixXd S(m, s);
	S.setZero();
	for (int i = 0; i < s; i++) {
		double eval = es.eigenvalues().real()[idx[i]];
		if (eval > 0) {
			S(i, i) = std::sqrt(eval);
		}
	}

	MatrixXd SI = S.transpose();
	for (int i = 0; i < s; i++) {
		SI(i, i) = 1 / SI(i, i);
	}

	MatrixXd U = A * V * SI;

	cout << "A: " << endl;
	cout << A << endl << endl;
	cout << "U: " << endl;
	cout << U << endl << endl;
	cout << "S: " << endl;
	cout << S << endl << endl;
	cout << "SI:" << endl;
	cout << SI << endl << endl;
	cout << "VT: " << endl;
	cout << V.transpose() << endl << endl;
	cout << "U * S * VT: " << endl;
	cout << U * S * V.transpose() << endl << endl;


	// Thin SVD
	MatrixXd U2 = U;
	MatrixXd S2 = S.block(0, 0, m, m);
	MatrixXd VT2 = V.transpose().block(0, 0, m, s);

	cout << "Thin SVD:" << endl;
	cout << U2 * S2 * VT2 << endl << endl;


	// Compact SVD
	for (int m2 = 5; m2 > 0; m2--) {
		MatrixXd U3 = U.block(0, 0, m, m2);
		MatrixXd S3 = S.block(0, 0, m2, m2);
		MatrixXd VT3 = V.transpose().block(0, 0, m2, s);
		MatrixXd Ap = U3 * S3 * VT3;
		to_unsigned_char_range(Ap);
		cout << "Compact SVD(rank " << m2 << "):" << endl;
		cout << Ap << endl << endl;
	}

	return 0;
}