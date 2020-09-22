//// Author: DY Suh
//// Date : June 23, 2020
////
////  Gradient descend for line fitting
////  input : data1.txt     random data(almost linear) --> not working
////  input : 4 samples of the blank constructor --> working
////
////  Optimize for (a, b) of y = ax + b
////
//#include <iostream> // for cout
//#include <iomanip> // for setw()
//using namespace std;
//#include "GDlinearFn.h"
//
//float EE(float x0, float y0, float x1, float y1)
//{
//	return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1));
//}
//
//int main() {
//	linearFn ab;
//
//	cout << "distance btw y = x and (1, 3): ";
//	cout << ab.distance(1, 1, 1, 3) << endl;
//
//	cout << "sum of distances btw y = x and ";
//	cout << "{(0.1, 1.1), (-0.1, 0.9), (1.1, 2.1), (0.9, 1.9)}: ";
//	cout << ab.LossFn(1, 1) << endl;
//
//	float psi = 0.01, eta = 0.0000001; // for 4 sample (blank constructor)
//	float da = 0.01, db = 0.01;
//	float a0 = -2, b0 = 2;
//	float a1 = 2.1, b1 = -0.8;  // answer (a, b) = (1, 1)   blank constructor
//	int iteration = 0;
//
//	while (EE(a0, b0, a1, b1) > eta && iteration < 1000) {
//		a0 = a1;
//		b0 = b1;
//		a1 -= psi * ab.dfabda(a0, b0, da);
//		b1 -= psi * ab.dfabdb(a0, b0, db);
//		iteration++;
//	}
//	cout << iteration << "-th  E = " << EE(a0, b0, a1, b1) << endl;
//	cout << "a1 = " << a1 << ", b1 = " << b1 << endl;
//
//	return 0;
//}
//
//
//
//
////float gaussian(float x, float y, float mux, float muy, float sigx, float sigy, float peak)
////{
////	return (peak * exp(-pow((x - mux) / sigx, 2.0) - pow((y - muy) / sigy, 2.0)));
////}
////
////float fxy(float x, float y)
////{
////	return (gaussian(x, y, 1., 1., 1., 2., 4) + gaussian(x, y, -1., -1., 1., 1., 2));
////}
//
////#include <iostream>
////#include <cmath>
////using namespace std;
////
////float PI = 3.141592;
////
////float fxy(float x, float y) { return sin(2 * PI * x) * sin(4 * PI * y); }
////float dfxydx(float x, float dx, float y) { return (fxy(x + dx, y) - fxy(x, y)) / dx; }
////float dfxydy(float x, float y, float dy) { return (fxy(x, y + dy) - fxy(x, y)) / dy; }
////float EE(float x0, float y0, float x1, float y1) { return sqrt((x0 - x1) * (x0 - x1) + (y0 - y1) * (y0 - y1)); }
////
////
////int main()
////{
////	// sprint1
////	for (float x = 0.25; x <= 0.75; x += 0.5) {
////		for (float y = 0.125; y <= 0.875; y += 0.25) {
////			cout << fxy(x, y) << endl;
////		}
////	}
////
////	float psi = 0.001, eta = 0.0005, dx = 0.01, dy = 0.01;
////	float x0 = 0.8, y0 = 0.8, x1 = 0.3, y1 = 0.3;
////
////	int iter = 0;
////	while (EE(x0, y0, x1, y1) > eta && iter < 100) {
////		x0 = x1;
////		y0 = y1;
////		x1 = x1 - psi * dfxydx(x0, dx, y0);
////		y1 = y1 - psi * dfxydy(x0, y0, dy);
////		iter++;
////	}
////
////	cout << iter << "-th " << x1 << ' ' << y1 << endl;
////
////	return 0;
////}