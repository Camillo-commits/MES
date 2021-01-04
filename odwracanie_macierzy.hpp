#ifndef MACIERZ
#define MACIERZ
#include <iostream>
#include <string>
#include <conio.h>
#include <fstream>
#include <cmath>

void pivoting(float ** arr,float ** arr2, int i, int number_of_rows, int number_of_columns);

void macierz_odwr(float** arr, int number_of_columns,int number_of_rows){
	std::fstream file;
	std::string text;

	int number_of_x = number_of_columns;
	float ** arr2;
	
	float float_tmp;

	


		//sprawdzanie poprawnosci danych(nie sprawdzam liniowej niezaleznosci)
		
		//generowanie tablicy
		arr2 = new float *[number_of_rows];
		for (int i = 0; i < number_of_rows; ++i) {
			arr2[i] = new float[number_of_x /*+ 1*/];
		}

		//generowanie macierzy diagonalnej
		for (int i = 0; i < number_of_rows; ++i) {
			for (int j = 0; j < number_of_columns; ++j) {
				arr2[i][j] = 0;
			}
		}
		for (int i = 0; i < number_of_columns; ++i) {
			arr2[i][i] = 1;
		}
		
		//weryfikacja poprawnosci danych
	
		
		//diagonalizacja macierzy
		int iteration = number_of_rows - 1;
		for (int i = 0; i < number_of_rows; ++i) {
			//zamiana kolumn

			pivoting(arr,arr2, i, number_of_rows, number_of_x);

			//odejmowanie
			for (int k = 0; k < iteration; ++k) {
				//obliczenie mnoznika

				float mnoznik = (arr[k + i + 1][i] / arr[i][i]);
				
				//odjecie poszczegolnych kolumn
				for (int j = 0; j < number_of_x /*+ 1*/; ++j) {
					arr[k + i + 1][j] -= arr[i][j] * mnoznik;
					arr2[k + i + 1][j] -= arr2[i][j] * mnoznik;
				}

			}
			--iteration;
			
		}
		//po dolnej
	
		//usuniecie czesci gornej macierzy trojkatnej
		iteration = number_of_rows;//- 1;
		for (int i = number_of_columns; i > 0;  --i) {
			for (int j = 1; j < iteration; ++j) {
				//		  w tej samej kol jeden wy¿ej   //ta na diagonali
				float mnoznik = (arr[i - 1 - j][i - 1] / arr[i - 1][i - 1]);
				arr[i - 1 - j][i - 1] -= arr[i - 1][i - 1] * mnoznik;
				arr2[i - 1 - j][i - 1] -= arr2[i - 1][i - 1] * mnoznik;
				
			}
			--iteration;
		}
	
		//doprowadzenie wartosci na diagonali do 1
		for (int i = 0; i < number_of_columns; ++i) {
			arr2[i][i] *= (1 / arr[i][i]);
			arr[i][i] *= (1 / arr[i][i]);
		}
		//przepisanie arr2 do arr
		
		for (int i = 0; i < number_of_rows; ++i) {
			for (int j = 0; j < number_of_columns; ++j) {
				arr[i][j] = arr2[i][j];
			}
		}

		

		//obliczenie kazdej ze zmiennych 
		/*
		float *result_arr = new float[number_of_x];
		//ostatnia zmienna
		for (int i = 0; i < number_of_x; ++i) {
			result_arr[i] = 0;
		}
		result_arr[number_of_x - 1] = (arr[number_of_rows - 1][number_of_x] / arr[number_of_rows - 1][number_of_x - 1]);



		for (int i = number_of_rows - 1; i > 0; --i) {
			float suma = 0;
			//suma wszystkich zmiennych pomnozonych o ich wsp bez szukanej
			for (int j = 0; j < number_of_x; ++j) {
				if ((i - 1) != j) {
					suma += (arr[i - 1][j] * result_arr[j]);
				}
			}
			result_arr[i - 1] = ((arr[i - 1][number_of_x] - suma) / arr[i - 1][i - 1]);
		}
		//podanie wyniku
		std::cout << "Wynik: " << std::endl;
		for (int i = 0; i < number_of_x; ++i) {
			std::cout << "X" << i << " = " << result_arr[i] << std::endl;
		}*/

	
	return;
}


void pivoting(float ** arr,float** arr2, int i, int number_of_rows, int number_of_columns) {
	//i pokazuje od ktorego wiersza zaczac i od jakiej kolumny
	int x = i;
	int y = i;
	int max = fabs(arr[y][x]);
	int max_row = y;
	float* tmp_arr = new float[number_of_columns /*+ 1*/];
	float* tmp_arr2 = new float[number_of_columns /*+ 1*/];
	//znalezienie maxa
	for (int j = 0; j < number_of_rows - i; ++j) {
		if (fabs(arr[y][x]) > max) {
			max = fabs(arr[y][x]);
			max_row = y;
		}
		++y;
	}
	//zapisanie zmienianej kolumny
	for (int j = 0; j < number_of_columns /*+ 1*/; ++j) {
		tmp_arr[j] = arr[i][j];
		tmp_arr2[j] = arr2[i][j];
	}
	//nadpisanie kolumn
	for (int j = 0; j < number_of_columns /*+ 1*/; ++j) {
		arr[i][j] = arr[max_row][j];
		arr2[i][j] = arr2[max_row][j];
	}
	for (int j = 0; j < number_of_columns /*+ 1*/; ++j) {
		arr[max_row][j] = tmp_arr[j];
		arr2[max_row][j] = tmp_arr2[j];
	}
	//
	delete[] tmp_arr;
	delete[] tmp_arr2;
	
}
#endif // !MACIERZ
