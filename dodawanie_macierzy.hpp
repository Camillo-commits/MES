#ifndef DODAWANIE_MACIERZY
#define DODAWANIE_MACIERZY
#include <iostream>

float** dodanie_macierzy(float** A, float** B, int ilosc_kolumn, int ilosc_wierszy) {
	float** C = new float*[ilosc_wierszy];
	for (int i = 0; i < ilosc_wierszy; ++i) {
		C[i] = new float[ilosc_kolumn];
	}
	//zerowanie macierzy C
	/*for (int i = 0; i < ilosc_wierszy; ++i) {
		for (int j = 0; j < ilosc_kolumn; ++j) {
			C[i][j] = 0;
		}
	}*/

	for (int i = 0; i < ilosc_wierszy; ++i) {
		for (int j = 0; j < ilosc_kolumn; ++j) {
			C[i][j] = (A[i][j] + B[i][j]);
		}
	}
	return C;
}



#endif // !DODAWANIE_MACIERZY
