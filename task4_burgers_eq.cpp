// heat_lab_5.cpp : ???? ???? ???????? ??????? "main". ????? ?????????? ? ????????????? ?????????? ?????????.
//


#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
using namespace std;

int main() {
	int n = 101, iter = 0;
	double dx = 1.0 / (n - 1),dy = 1.0 / (n - 1), dt, dif = 0.0, eps = pow(10,-6);
	double **u0 = new double *[n];
	double **u = new double *[n];
	for (int i = 0; i < n; i++) {
		u0[i] = new double[n];
		u[i] = new double[n];
	}
	dt = 0.5 / (1.0 / (dx * dx) + 1.0 / (dy * dy));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			u0[i][j] = 0.0;
			u[i][j] = 0.0;
		}
	}

	for (int i = 0; i < n; i++)
	{
		u0[i][0] = 1.0;
		u[i][0] = 1.0;

	}
		do {
		for (int i = 1; i < n - 1; i++)
			for (int j = 1; j < n - 1; j++)
				u[i][j] = u0[i][j] + dt * ((u0[i + 1][j] - 2. * u0[i][j] + u0[i - 1][j]) / (dx * dx) +
				(u0[i][j + 1] - 2. * u0[i][j] + u0[i][j - 1]) / (dy * dy));

		dif = 0.0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				if (dif< abs(u[i][j] - u0[i][j])) {
					dif = abs(u[i][j] - u0[i][j]);
				}
			}
		}
	

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				u0[i][j] = u[i][j];
			}
		}

		iter++;
	} while (dif > eps);


	ofstream fout("task5.dat");
fout << "VARIABLES = \"X\",\"Y\",\"u\"" << endl;
	fout << "ZONE I=" << n << ",J=" << n<< ",F=POINT" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			fout << i * dx << '\t' << j * dy << '\t' << u[i][j] << endl;
		}
	}

	cout << "max difference: " << dif << endl;
	cout << "number of iterations: " << iter << endl;

	return 0;
}
