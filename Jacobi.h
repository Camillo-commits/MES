#ifndef JACOBI
#define JACOBI

#include "odwracanie_macierzy.hpp"
#include "mnozenie_macierzy.hpp"
#include "dodawanie_macierzy.hpp"
#include<iostream>

void changeRows(double** arr, int row1, int row2, int sizeX);
void changeRows(double* arr, int row1, int row2);

double* resolveSOEJacobi(double** params, double* equals,int size) {
	double** X = new double*[size];
	double* B = new double[size];
	double* alfa = new double[size];

	for (int i = 0; i < size; ++i) {
		X[i] = new double[size];
	}

	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			X[i][j] = params [i][j];
		}
		B[i] = equals[i];
	}

	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			if (X[j][j] == 0) {
				bool ok = false;
				for (int k = 0; k < size; ++k) {
					if (X[k][j] != 0) {
						changeRows(X, j, k, size);
						changeRows(B, j, k);
						ok = true;
					}
				}
				if (!ok) {
					throw "Not solvable";
				}
			}
		}
		//calcAlfa
		for (int j = 0; j < size; ++j) {
			if (i != j)
				alfa[j] = X[j][i] / X[i][i];
			else
				alfa[j] = 0;
		}
		
		//calculate
		
		for (int j = 0; j < size; ++j) {
			for (int k = 0; k < size; ++k) {
				X[j][k] -= alfa[j] * X[i][k];
			}
			B[j] -= alfa[j] * B[i];
		}
		/*for (int j = 0; j < size; ++j) {
			for (int k = 0; k < size; ++k) {
				std::cout << X[j][k] << "		";
			}
			std::cout << " | " << B[j] << std::endl;
		}*/
	}
	for (int i = 0; i < size; ++i) {
		if (X[i][i] != 1) {
			B[i] /= X[i][i];
			X[i][i] = 1;
		}
	}
	return B;

}


void changeRows(double* arr, int row1, int row2) {
	double tmp;
	
	tmp = arr[row1];
	arr[row1] = arr[row2];
	arr[row2] = tmp;
}

void changeRows(double** arr, int row1, int row2, int sizeX) {
	double tmp;
	for (int i = 0; i < sizeX; ++i) {
		tmp = arr[row1][i];
		arr[row1][i] = arr[row2][i];
		arr[row2][i] = tmp;
	}
}

#endif // !JACOBI
