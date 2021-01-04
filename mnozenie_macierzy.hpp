#ifndef MNOZENIE
#define MNOZENIE
#include <iostream>

float** pomnozA_razy_B(float** A, float** B, int ilosc_wierszy_A,int ilosc_kolumn_A, int ilosc_wierszy_B, int ilosc_kolumn_B) {
	//indeksowanie najpierw wspolrzedna Y potem X
	if (ilosc_kolumn_A != ilosc_wierszy_B) {
		std::cout << "mnozenie nie mozliwe" << std::endl;
		
	}
	//tworzenie macierzy wynikowej
	float** C = new float* [ilosc_wierszy_B];
	for (int i = 0; i < ilosc_wierszy_B; ++i) {
		C[i] = new float[ilosc_kolumn_A];
	}
	int j = 0;
	int k = 0;
	float float_tmp = 0;

	//wymnozenie jednego wiersza
	for (int y = 0; y < ilosc_wierszy_B; ++y) {
		for (int x = 0; x < ilosc_kolumn_A; ++x) {
			//obliczenie pojedynczej komorki C
			for (int l = 0, i = 0; l < ilosc_wierszy_B; ++l) {
				//std::cout << "macierz A : " << j << " " << i  << " wartosc: " << A[j][i] << std::endl;
				//std::cout << "macierz B : " << l << " " << k << " wartosc: " << B[l][k] << std::endl;

				float_tmp += A[j][i] * B[l][k];
				//std::cout << "float_tmp: " << float_tmp << std::endl << std::endl;
				++i;
			}
			C[j][k] = float_tmp;
			++k;
			float_tmp = 0;
		}
		k = 0;
		++j;
	}
	//std::cout << "wynik: " << std::endl;
	//for (int y = 0; y < ilosc_wierszy_B; ++y) {
	//	for (int x = 0; x < ilosc_kolumn_A; ++x) {
	//		std::cout << C[y][x] << "  ";
	//	}
	//	std::cout << std::endl;
	//}
	return C;
}

#endif // !MNOZENIE
