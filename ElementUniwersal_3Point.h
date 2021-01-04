#ifndef ELEM_UNIWERSAL_3_Point
#define ELEM_UNIWERSAL_3_Point

#include <iostream>

class ElemUniwersal4_3point {
private:
	double pierw3_5 = 0.7745966692414834;
	
	//								1 / sqrt(3) ( punkty calkowania )
	double E[9] = { -pierw3_5, 0, pierw3_5, -pierw3_5, 0, pierw3_5, -pierw3_5, 0, pierw3_5 }; //like X
	double n[9] = { -pierw3_5, -pierw3_5,  -pierw3_5, 0, 0, 0, pierw3_5, pierw3_5, pierw3_5}; //like Y
	//
	//jakobiany dla kazdego z 9ciu punktow
	double J1[2][2];
	double J2[2][2];
	double J3[2][2];
	double J4[2][2];
	double J5[2][2];
	double J6[2][2];
	double J7[2][2];
	double J8[2][2];
	double J9[2][2];

	//jakobiany odwrotne
	double revJ1[2][2];
	double revJ2[2][2];
	double revJ3[2][2];
	double revJ4[2][2];
	double revJ5[2][2];
	double revJ6[2][2];
	double revJ7[2][2];
	double revJ8[2][2];
	double revJ9[2][2];

	//wyznaczniki macierzy
	double detJ1;
	double detJ2;
	double detJ3;
	double detJ4;
	double detJ5;
	double detJ6;
	double detJ7;
	double detJ8;
	double detJ9;

	//dN/dE
	double derv_dN_dE[4][9];
	//dN/dn
	double derv_dN_dn[4][9];

	double dn_dxy_1[2][4];
	double dn_dxy_2[2][4];
	double dn_dxy_3[2][4];
	double dn_dxy_4[2][4];
	double dn_dxy_5[2][4];
	double dn_dxy_6[2][4];
	double dn_dxy_7[2][4];
	double dn_dxy_8[2][4];
	double dn_dxy_9[2][4];

	//wagi calkowania
	double w[3] = { 0.5555555555555555 , 0.8888888888888888, 0.5555555555555555 };

	double H[4][4];

	double k_t = 25;
	double alfa = 300;
	double t_inf = 1200;

	double C[4][4];
	double H_BC[4][4];
	double N[4][9];
	double P[4];

	void cTmpTimesDet(double** cTmp, int detNumber);
	void cTmpTimesW(double** cTmp, int NNumber);
	void addToC(double** cTmp);
	void H_BCTimesDet(int wall, double det);
	void H_BCTimesAlfa();

	void PTimesAlfaT_inf();
	void gen_dn_dxy(double revJ[2][2], double tmp[2][4], double dn_dxy[2][4]);
public:	
	void genJacobians(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
	void genC(double c,double ro);
	void jClear(int numOfJ);
	void jPrint(int numOfJ);
	void derv_dN_dEPrint();
	void derv_dN_dnPrint();
	void revJPrint(int numOfJ);
	void dn_xyPrint(int numOfdx_dx);
	void hPrint();
	void CClear();
	void CPrint();
	void genH_BC(double x1, double y1, int BC1, double x2, double y2, int BC2, double x3, double y3, int BC3, double x4, double y4, int BC4);
	void H_BCClear();
	void H_BCPrint();
	
	void genP(double x1, double y1, int BC1, double x2, double y2, int BC2, double x3, double y3, int BC3, double x4, double y4, int BC4);
	void PClear();
	void PPrint();
	double getAlfa() { return alfa; }
	double getT_inf() { return t_inf; }
	void setAlfa(double Alfa) { this->alfa = Alfa; }
	void setT_inf(double T_inf) { this->t_inf = T_inf; }
	double getK() { return k_t; }
	void setK(double K) { k_t = K; }

	void addH_BC2H();

	double** getH();
	double** getC();
	double* getP();
	ElemUniwersal4_3point() = default;
};


void ElemUniwersal4_3point::addH_BC2H() {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			H[i][j] += H_BC[i][j];
		}
	}
}

double* ElemUniwersal4_3point::getP(){
	double* tmp = new double[4];
	for (int i = 0; i < 4; ++i) {
		tmp[i] = P[i];
	}
	return tmp;
}

void ElemUniwersal4_3point::PClear() {
	for (int i = 0; i < 4; ++i) {
		P[i] = 0;
	}
}

void ElemUniwersal4_3point::PPrint() {
	std::cout << std::endl << "P: " << std::endl;
	for (int i = 0; i < 4; ++i) {
		std::cout << P[i] << "		";
	}
}

void ElemUniwersal4_3point::genP(double x1, double y1, int BC1, double x2, double y2, int BC2, double x3, double y3, int BC3, double x4, double y4, int BC4) {
	PClear();
	double N[3][2];
	double L;
	for (int i = 0; i < 3; ++i) {
		N[i][0] = (1 - E[i]) / 2;
		N[i][1] = (1 + E[i]) / 2;
	}

	//dolny bok
	if (BC1 == 1 && BC2 == 1) {
		L = (x2 - x1) / 2;
		for (int i = 0; i < 3; ++i) {
			P[0] += N[i][0] * L * w[i];
			P[1] += N[i][1] * L * w[i];
		}
		//PTimesDet(0, L);
	}
	//prawy bok
	if (BC2 == 1 && BC3 == 1) {
		L = (y3 - y2) / 2;
		for (int i = 0; i < 3; ++i) {
			P[1] += N[i][0] * L * w[i];
			P[2] += N[i][1] * L * w[i];
		}
		//PTimesDet(1, L);
	}
	//gorny bok
	if (BC3 == 1 && BC4 == 1) {
		L = (x3 - x4) / 2;
		for (int i = 0; i < 3; ++i) {
			P[2] += N[i][0] * L * w[i];
			P[3] += N[i][1] * L * w[i];
		}
		//PTimesDet(2, L);
	}
	//lewy bok
	if (BC4 == 1 && BC1 == 1) {
		L = (y4 - y1) / 2;
		for (int i = 0; i < 3; ++i) {
			P[3] += N[i][0] * L * w[i];
			P[0] += N[i][1] * L * w[i];
		}
		//PTimesDet(3, L);
	}

	PTimesAlfaT_inf();
}

void ElemUniwersal4_3point::PTimesAlfaT_inf() {
	for (int i = 0; i < 4; ++i) {
		P[i] *= -alfa;
		P[i] *= t_inf;
	}
}

void ElemUniwersal4_3point::genH_BC(double x1, double y1, int BC1, double x2, double y2, int BC2, double x3, double y3, int BC3, double x4, double y4, int BC4) {
	H_BCClear();
	double N[3][2];
	double L;
	for (int i = 0; i < 3; ++i) {
		N[i][0] = (1 - E[i]) / 2;
		N[i][1] = (1 + E[i]) / 2;
	}	
	double tmp1[2][2], tmp2[2][2], tmp3[2][2], tmp4[2][2];
	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < 2; ++j) {
			tmp1[i][j] = 0;
			tmp2[i][j] = 0;
			tmp3[i][j] = 0;
			tmp4[i][j] = 0;
		}
	}
	//dolny bok
	if (BC1 == 1 && BC2 == 1) {
		L = (x2 - x1) / 2;
		for (int i = 0; i < 3; ++i) {
			tmp1[0][0] += N[i][0] * N[i][0] * w[i];
			tmp1[0][1] += N[i][0] * N[i][1] * w[i];
			tmp1[1][0] += N[i][0] * N[i][1] * w[i];
			tmp1[1][1] += N[i][1] * N[i][1] * w[i];
			
		}
		tmp1[0][0] *= L;
		tmp1[0][1] *= L;
		tmp1[1][0] *= L;
		tmp1[1][1] *= L;

		//H_BCTimesDet(0, L);
	}
	//prawy bok
	if (BC2 == 1 && BC3 == 1) {
		L = (y3 - y2) / 2;
		for (int i = 0; i < 3; ++i) {
			tmp2[0][0] += N[i][0] * N[i][0] * w[i];
			tmp2[0][1] += N[i][0] * N[i][1] * w[i];
			tmp2[1][0] += N[i][0] * N[i][1] * w[i];
			tmp2[1][1] += N[i][1] * N[i][1] * w[i];
		}
		//H_BCTimesDet(1, L);
		tmp2[0][0] *= L;
		tmp2[0][1] *= L;
		tmp2[1][0] *= L;
		tmp2[1][1] *= L;
	}
	//gorny bok
	if (BC3 == 1 && BC4 == 1) {
		L = (x3 - x4) / 2;
		for (int i = 0; i < 3; ++i) {
			tmp3[0][0] += N[i][0] * N[i][0] * w[i];
			tmp3[0][1] += N[i][0] * N[i][1] * w[i];
			tmp3[1][0] += N[i][0] * N[i][1] * w[i];
			tmp3[1][1] += N[i][1] * N[i][1] * w[i];
		}
		//H_BCTimesDet(2, L);
		tmp3[0][0] *= L;
		tmp3[0][1] *= L;
		tmp3[1][0] *= L;
		tmp3[1][1] *= L;
	}
	//lewy bok
	if (BC4 == 1 && BC1 == 1) {
		L = (y4 - y1) / 2;
		for (int i = 0; i < 3; ++i) {
			tmp4[0][0] += N[i][0] * N[i][0] * w[i];
			tmp4[0][1] += N[i][0] * N[i][1] * w[i];
			tmp4[1][0] += N[i][0] * N[i][1] * w[i];
			tmp4[1][1] += N[i][1] * N[i][1] * w[i];
		}
		//H_BCTimesDet(3, L);
		tmp4[0][0] *= L;
		tmp4[0][1] *= L;
		tmp4[1][0] *= L;
		tmp4[1][1] *= L;
	}

	H_BC[0][0] += tmp1[0][0];
	H_BC[0][1] += tmp1[0][1];
	H_BC[1][0] += tmp1[1][0];
	H_BC[1][1] += tmp1[1][1];

	H_BC[1][1] += tmp2[0][0];
	H_BC[1][2] += tmp2[0][1];
	H_BC[2][1] += tmp2[1][0];
	H_BC[2][2] += tmp2[1][1];

	H_BC[2][2] += tmp3[0][0];
	H_BC[2][3] += tmp3[0][1];
	H_BC[3][2] += tmp3[1][0];
	H_BC[3][3] += tmp3[1][1];

	H_BC[0][0] += tmp4[0][0];
	H_BC[0][3] += tmp4[0][1];
	H_BC[3][0] += tmp4[1][0];
	H_BC[3][3] += tmp4[1][1];

	H_BCTimesAlfa();
}

void ElemUniwersal4_3point::H_BCPrint() {
	std::cout << std::endl << "H_BC" << std::endl;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << H_BC[i][j] << "		";
		}
		std::cout << std::endl;
	}
}

void ElemUniwersal4_3point::H_BCClear() {

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			H_BC[i][j] = 0;
		}
	}

}

void ElemUniwersal4_3point::H_BCTimesDet(int wall, double det) {
	if (wall == 0) {
		H_BC[0][0] *= det;
		H_BC[0][1] *= det;
		H_BC[1][0] *= det;
		H_BC[1][1] *= det;
	}
	if (wall == 1) {
		H_BC[1][1] *= det;
		H_BC[1][2] *= det;
		H_BC[2][1] *= det;
		H_BC[2][2] *= det;
	}
	if (wall == 2) {
		H_BC[2][2] *= det;
		H_BC[2][3] *= det;
		H_BC[3][2] *= det;
		H_BC[3][3] *= det;
	}
	if (wall == 3) {
		H_BC[0][0] *= det;
		H_BC[0][3] *= det;
		H_BC[3][0] *= det;
		H_BC[3][3] *= det;
	}
}

void ElemUniwersal4_3point::H_BCTimesAlfa() {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			H_BC[i][j] *= alfa;
		}
	}
}

void ElemUniwersal4_3point::CPrint() {
	std::cout << std::endl << "C:" << std::endl;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << C[i][j] << "			";
		}
		std::cout << std::endl;
	}
}

void ElemUniwersal4_3point::CClear() {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			C[i][j] = 0;
		}
	}
}

void ElemUniwersal4_3point::addToC(double** cTmp) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			C[i][j] += cTmp[i][j];
		}
	}
}

void ElemUniwersal4_3point::cTmpTimesW(double** cTmp, int NNumber) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			switch (NNumber) {
			case 0:
				cTmp[i][j] *= w[0] * w[0] ;
				break;
			case 1:
				cTmp[i][j] *= w[0] * w[1] ;
				break;
			case 2:
				cTmp[i][j] *= w[0] * w[2] ;
				break;
			case 3:
				cTmp[i][j] *= w[1] * w[0] ;
				break;
			case 4:
				cTmp[i][j] *= w[1] * w[1] ;
				break;
			case 5:
				cTmp[i][j] *= w[1] * w[2] ;
				break;
			case 6:
				cTmp[i][j] *= w[2] * w[0] ;
				break;
			case 7:
				cTmp[i][j] *= w[2] * w[1] ;
				break;
			case 8:
				cTmp[i][j] *= w[2] * w[2] ;
				break;
			default:
				break;
			}
		}
	}
}


void ElemUniwersal4_3point::cTmpTimesDet(double** cTmp, int detNumber) {
	for (int i = 0; i < 4; ++i) {
		for(int j = 0;j<4;++j){
			switch (detNumber) {
			case 0:
				cTmp[i][j] *= detJ1;
				break;
			case 1:
				cTmp[i][j] *= detJ2;
				break;
			case 2:
				cTmp[i][j] *= detJ3;
				break;
			case 3:
				cTmp[i][j] *= detJ4;
				break;
			case 4:
				cTmp[i][j] *= detJ5;
				break;
			case 5:
				cTmp[i][j] *= detJ6;
				break;
			case 6:
				cTmp[i][j] *= detJ7;
				break;
			case 7:
				cTmp[i][j] *= detJ8;
				break;
			case 8:
				cTmp[i][j] *= detJ9;
				break;
			default:
				break;
			}
		}
	}
}

void ElemUniwersal4_3point::genC(double c, double ro) {
	
	CClear();

	for (int i = 0; i < 9; ++i) {
		N[0][i] = 0.25 * (1 - E[i]) * (1 - n[i]);
		N[1][i] = 0.25 * (1 + E[i]) * (1 - n[i]);
		N[2][i] = 0.25 * (1 + E[i]) * (1 + n[i]);
		N[3][i] = 0.25 * (1 - E[i]) * (1 + n[i]);
	}
	/*
	for (int i = 0; i < 9; ++i) {
		std::cout << i << std::endl;
		std::cout << N[0][i] << std::endl;
		std::cout << N[1][i] << std::endl;
		std::cout << N[2][i] << std::endl;
		std::cout << N[3][i] << std::endl;
	}
	*/
	double** cTmp = new double*[4];
	for (int i = 0; i < 4; ++i) {
		cTmp[i] = new double[4];
	}

	
	//for each row
	for (int i = 0; i < 9; ++i) {
		//clear cTmp
		for (int k = 0; k < 4; ++k) {
			for (int j = 0; j < 4; ++j) {
				cTmp[j][k] = 0;
			}
		}

		for (int a = 0; a < 4; ++a) {
			for (int b = 0; b < 4; ++b) {
				cTmp[b][a] = N[a][i] * N[b][i];// * c * ro;
				double tmp = cTmp[b][a];
			}
		}
		/*
		std::cout << i << "wart c" << std::endl;
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << cTmp[i][j] << "		";
			}
			std::cout << std::endl;
		}*/
		cTmpTimesDet(cTmp, i);
		/*std::cout << i << "wart c po det" << std::endl;
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << cTmp[i][j] << "		";
			}
			std::cout << std::endl;
		}*/
		cTmpTimesW(cTmp, i);
		/*std::cout << i << "wart c po W" << std::endl;
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << cTmp[i][j] << "		";
			}
			std::cout << std::endl;
		}*/
		addToC(cTmp);
		
		
		//std::cout << std::endl;
		
	}
	
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			C[i][j] *= c * ro;// * detJ1;
		}
	}
	
}



double** ElemUniwersal4_3point::getH() {
	double **tmp = new double*[4];
	for (int i = 0; i < 4; ++i) {
		tmp[i] = new double[4];
	}
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			tmp[j][i] = H[j][i];
		}
	}
	return tmp;
}


double** ElemUniwersal4_3point::getC() {
	double **tmp = new double*[4];
	for (int i = 0; i < 4; ++i) {
		tmp[i] = new double[4];
	}
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			tmp[j][i] = C[j][i];
		}
	}
	return tmp;
}

void ElemUniwersal4_3point::hPrint() {
	std::cout << std::endl << "H" << std::endl;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << H[j][i] << " ";
		}
		std::cout << std::endl;
	}

}


void ElemUniwersal4_3point::revJPrint(int numOfJ) {
	std::cout << std::endl << "revJ" << numOfJ << std::endl;
	switch (numOfJ) {

	case 1:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				std::cout << revJ1[j][i] << " ";
			}
			std::cout << std::endl;
		}
		break;
	case 2:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				std::cout << revJ2[j][i] << " ";
			}
			std::cout << std::endl;
		}
		break;
	case 3:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				std::cout << revJ3[j][i] << " ";
			}
			std::cout << std::endl;
		}
		break;
	case 4:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				std::cout << revJ4[j][i] << " ";
			}
			std::cout << std::endl;
		}
	case 5:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				std::cout << revJ5[j][i] << " ";
			}
			std::cout << std::endl;
		}
		break;
	case 6:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				std::cout << revJ6[j][i] << " ";
			}
			std::cout << std::endl;
		}
		break;
	case 7:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				std::cout << revJ7[j][i] << " ";
			}
			std::cout << std::endl;
		}
		break;
	case 8:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				std::cout << revJ8[j][i] << " ";
			}
			std::cout << std::endl;
		}
		break;
	case 9:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 2; ++j) {
				std::cout << revJ9[j][i] << " ";
			}
			std::cout << std::endl;
		}
		break;
	}
}

void ElemUniwersal4_3point::dn_xyPrint(int numOfdx_dy) {
	std::cout << std::endl << "dn_xy" << numOfdx_dy << std::endl;
	switch (numOfdx_dy) {

	case 1:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << dn_dxy_1[i][j] << " ";
			}
			std::cout << std::endl;
		}
		break;
	case 2:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << dn_dxy_2[i][j] << " ";
			}
			std::cout << std::endl;
		}
		break;
	case 3:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << dn_dxy_3[i][j] << " ";
			}
			std::cout << std::endl;
		}
		break;
	case 4:
		for (int i = 0; i < 2; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << dn_dxy_4[i][j] << " ";
			}
			std::cout << std::endl;
		}
		break;
	default:
		break;
	}
}

void ElemUniwersal4_3point::jPrint(int numOfJ) {
	std::cout << std::endl << "J" << numOfJ << std::endl;
	switch (numOfJ)
	{
	case 1:
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J1[i][0] << "      ";
		}
		std::cout << std::endl;
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J1[i][1] << "      ";
		}
		break;
	case 2:
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J2[i][0] << "      ";
		}
		std::cout << std::endl;
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J2[i][1] << "      ";
		}
		break;
	case 3:
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J3[i][0] << "      ";
		}
		std::cout << std::endl;
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J3[i][1] << "      ";
		}
		break;
	case 4:
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J4[i][0] << "      ";
		}
		std::cout << std::endl;
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J4[i][1] << "      ";
		}
		break;
	case 5:
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J5[i][0] << "      ";
		}
		std::cout << std::endl;
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J5[i][1] << "      ";
		}
		break;

	case 6:
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J6[i][0] << "      ";
		}
		std::cout << std::endl;
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J6[i][1] << "      ";
		}
		break;
	case 7:
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J7[i][0] << "      ";
		}
		std::cout << std::endl;
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J7[i][1] << "      ";
		}
		break;

	case 8:
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J8[i][0] << "      ";
		}
		std::cout << std::endl;
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J8[i][1] << "      ";
		}
		break;
	case 9:
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J9[i][0] << "      ";
		}
		std::cout << std::endl;
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J9[i][1] << "      ";
		}
		break;
	default:
		break;
	}
}

void ElemUniwersal4_3point::jClear(int numOfJ) {
	switch (numOfJ)
	{
	case 1:
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				J1[j][i] = 0;
			}
		}
		break;
	case 2:
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				J2[j][i] = 0;
			}
		}
		break;
	case 3:
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				J3[j][i] = 0;
			}
		}
		break;
	case 4:
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				J4[j][i] = 0;
			}
		}
		break;
	case 5:
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				J5[j][i] = 0;
			}
		}
		break;
	case 6:
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				J6[j][i] = 0;
			}
		}
		break;
	case 7:
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				J7[j][i] = 0;
			}
		}
		break;
	case 8:
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				J8[j][i] = 0;
			}
		}
		break;
	case 9:
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				J9[j][i] = 0;
			}
		}
		break;
	default:
		break;
	}
}

void ElemUniwersal4_3point::derv_dN_dEPrint() {
	std::cout << std::endl;
	std::cout << "derv_dN/dE" << std::endl;
	for (int i = 0; i < 9; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << derv_dN_dE[j][i] << " ";
		}
		std::cout << std::endl;
	}


}

void ElemUniwersal4_3point::derv_dN_dnPrint() {
	std::cout << std::endl;
	std::cout << "derv_dN/dn" << std::endl;
	for (int i = 0; i < 9; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << derv_dN_dn[j][i] << " ";
		}
		std::cout << std::endl;
	}
}

void ElemUniwersal4_3point::gen_dn_dxy(double revJ[2][2], double tmp[2][4], double dn_dxy[2][4]) {
	int j = 0;
	int k = 0;
	float float_tmp = 0;

	//wymnozenie jednego wiersza
	for (int y = 0; y < 2; ++y) {
		for (int x = 0; x < 4; ++x) {
			//obliczenie pojedynczej komorki C
			for (int l = 0, i = 0; l < 2; ++l) {
				//std::cout << "macierz A : " << j << " " << i  << " wartosc: " << A[j][i] << std::endl;
				//std::cout << "macierz B : " << l << " " << k << " wartosc: " << B[l][k] << std::endl;

				float_tmp += revJ[j][i] * tmp[l][k];
				//std::cout << "float_tmp: " << float_tmp << std::endl << std::endl;
				++i;
			}
			dn_dxy[j][k] = float_tmp;
			++k;
			float_tmp = 0;
		}
		k = 0;
		++j;
	}
}

void ElemUniwersal4_3point::genJacobians(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
	double xArr[4] = { x1,x2,x3,x4 };
	double yArr[4] = { y1,y2,y3,y4 };

	for (int i = 0; i < 9; ++i) {
		//dN/dE
		derv_dN_dE[0][i] = -0.25 * (1 - n[i]);
		derv_dN_dE[1][i] = 0.25 * (1 - n[i]);
		derv_dN_dE[2][i] = 0.25 * (1 + n[i]);
		derv_dN_dE[3][i] = -0.25 * (1 + n[i]);
		//dN/dn
		derv_dN_dn[0][i] = -0.25 * (1 - E[i]);
		derv_dN_dn[1][i] = -0.25 * (1 + E[i]);
		derv_dN_dn[2][i] = 0.25 * (1 + E[i]);
		derv_dN_dn[3][i] = 0.25 * (1 - E[i]);
	}
	//J1
	jClear(1);
	jClear(2);
	jClear(3);
	jClear(4);
	jClear(5);
	jClear(6);
	jClear(7);
	jClear(8);
	jClear(9);
	for (int i = 0; i < 4; ++i) {
		//	dX/dE
		J1[0][0] += derv_dN_dE[i][0] * xArr[i];
		//	dY/dE
		J1[1][0] += derv_dN_dE[i][0] * yArr[i];
		//	dX/dn
		J1[0][1] += derv_dN_dn[i][0] * xArr[i];
		//	dY/dn
		J1[1][1] += derv_dN_dn[i][0] * yArr[i];
	
		//	dX/dE
		J2[0][0] += derv_dN_dE[i][1] * xArr[i];
		//	dY/dE
		J2[1][0] += derv_dN_dE[i][1] * yArr[i];
		//	dX/dn
		J2[0][1] += derv_dN_dn[i][1] * xArr[i];
		//	dY/dn
		J2[1][1] += derv_dN_dn[i][1] * yArr[i];
	
		//	dX/dE
		J3[0][0] += derv_dN_dE[i][2] * xArr[i];
		//	dY/dE
		J3[1][0] += derv_dN_dE[i][2] * yArr[i];
		//	dX/dn
		J3[0][1] += derv_dN_dn[i][2] * xArr[i];
		//	dY/dn
		J3[1][1] += derv_dN_dn[i][2] * yArr[i];
	
		//	dX/dE
		J4[0][0] += derv_dN_dE[i][3] * xArr[i];
		//	dY/dE
		J4[1][0] += derv_dN_dE[i][3] * yArr[i];
		//	dX/dn
		J4[0][1] += derv_dN_dn[i][3] * xArr[i];
		//	dY/dn
		J4[1][1] += derv_dN_dn[i][3] * yArr[i];

		//	dX/dE
		J5[0][0] += derv_dN_dE[i][4] * xArr[i];
		//	dY/dE
		J5[1][0] += derv_dN_dE[i][4] * yArr[i];
		//	dX/dn
		J5[0][1] += derv_dN_dn[i][4] * xArr[i];
		//	dY/dn
		J5[1][1] += derv_dN_dn[i][4] * yArr[i];
	
		//	dX/dE
		J6[0][0] += derv_dN_dE[i][5] * xArr[i];
		//	dY/dE
		J6[1][0] += derv_dN_dE[i][5] * yArr[i];
		//	dX/dn
		J6[0][1] += derv_dN_dn[i][5] * xArr[i];
		//	dY/dn
		J6[1][1] += derv_dN_dn[i][5] * yArr[i];

		//	dX/dE
		J7[0][0] += derv_dN_dE[i][6] * xArr[i];
		//	dY/dE
		J7[1][0] += derv_dN_dE[i][6] * yArr[i];
		//	dX/dn
		J7[0][1] += derv_dN_dn[i][6] * xArr[i];
		//	dY/dn
		J7[1][1] += derv_dN_dn[i][6] * yArr[i];

		//	dX/dE
		J8[0][0] += derv_dN_dE[i][7] * xArr[i];
		//	dY/dE
		J8[1][0] += derv_dN_dE[i][7] * yArr[i];
		//	dX/dn
		J8[0][1] += derv_dN_dn[i][7] * xArr[i];
		//	dY/dn
		J8[1][1] += derv_dN_dn[i][7] * yArr[i];
	
		//	dX/dE
		J9[0][0] += derv_dN_dE[i][8] * xArr[i];
		//	dY/dE
		J9[1][0] += derv_dN_dE[i][8] * yArr[i];
		//	dX/dn
		J9[0][1] += derv_dN_dn[i][8] * xArr[i];
		//	dY/dn
		J9[1][1] += derv_dN_dn[i][8] * yArr[i];
	}
	/////
	//revJ1
	detJ1 = (J1[0][0] * J1[1][1]) - (J1[0][1] * J1[1][0]);
	revJ1[0][0] = J1[1][1] * (1 / detJ1);
	revJ1[1][0] = -J1[1][0] * (1 / detJ1);
	revJ1[0][1] = -J1[0][1] * (1 / detJ1);
	revJ1[1][1] = J1[0][0] * (1 / detJ1);
	////

	//revJ2
	detJ2 = (J2[0][0] * J2[1][1]) - (J2[0][1] * J2[1][0]);
	revJ2[0][0] = J2[1][1] * (1 / detJ2);
	revJ2[1][0] = -J2[1][0] * (1 / detJ2);
	revJ2[0][1] = -J2[0][1] * (1 / detJ2);
	revJ2[1][1] = J2[0][0] * (1 / detJ2);
	////

	//revJ3
	detJ3 = (J3[0][0] * J3[1][1]) - (J3[0][1] * J3[1][0]);
	revJ3[0][0] = J3[1][1] * (1 / detJ3);
	revJ3[1][0] = -J3[1][0] * (1 / detJ3);
	revJ3[0][1] = -J3[0][1] * (1 / detJ3);
	revJ3[1][1] = J3[0][0] * (1 / detJ3);
	////

	//revJ4
	detJ4 = (J4[0][0] * J4[1][1]) - (J4[0][1] * J4[1][0]);
	revJ4[0][0] = J4[1][1] * (1 / detJ4);
	revJ4[1][0] = -J4[1][0] * (1 / detJ4);
	revJ4[0][1] = -J4[0][1] * (1 / detJ4);
	revJ4[1][1] = J4[0][0] * (1 / detJ4);
	/////
	/////
	//revJ5
	detJ5 = (J5[0][0] * J5[1][1]) - (J5[0][1] * J5[1][0]);
	revJ5[0][0] = J5[1][1] * (1 / detJ5);
	revJ5[1][0] = -J5[1][0] * (1 / detJ5);
	revJ5[0][1] = -J5[0][1] * (1 / detJ5);
	revJ5[1][1] = J5[0][0] * (1 / detJ5);
	////

	//revJ6
	detJ6 = (J6[0][0] * J6[1][1]) - (J6[0][1] * J6[1][0]);
	revJ6[0][0] = J6[1][1] * (1 / detJ6);
	revJ6[1][0] = -J6[1][0] * (1 / detJ6);
	revJ6[0][1] = -J6[0][1] * (1 / detJ6);
	revJ6[1][1] = J6[0][0] * (1 / detJ6);
	////

	//revJ7
	detJ7 = (J7[0][0] * J7[1][1]) - (J7[0][1] * J7[1][0]);
	revJ7[0][0] = J7[1][1] * (1 / detJ7);
	revJ7[1][0] = -J7[1][0] * (1 / detJ7);
	revJ7[0][1] = -J7[0][1] * (1 / detJ7);
	revJ7[1][1] = J7[0][0] * (1 / detJ7);
	////

	//revJ8
	detJ8 = (J8[0][0] * J8[1][1]) - (J8[0][1] * J8[1][0]);
	revJ8[0][0] = J8[1][1] * (1 / detJ8);
	revJ8[1][0] = -J8[1][0] * (1 / detJ8);
	revJ8[0][1] = -J8[0][1] * (1 / detJ8);
	revJ8[1][1] = J8[0][0] * (1 / detJ8);
	/////

	//ref J9
	detJ9 = (J9[0][0] * J9[1][1]) - (J9[0][1] * J9[1][0]);
	revJ9[0][0] = J9[1][1] * (1 / detJ9);
	revJ9[1][0] = -J9[1][0] * (1 / detJ9);
	revJ9[0][1] = -J9[0][1] * (1 / detJ9);
	revJ9[1][1] = J9[0][0] * (1 / detJ9);
	/////
	//dn_dxy_1
	double tmp[2][4] = {
		{derv_dN_dE[0][0],derv_dN_dE[1][0],derv_dN_dE[2][0],derv_dN_dE[3][0]},
		{derv_dN_dn[0][0], derv_dN_dn[1][0],derv_dN_dn[2][0],derv_dN_dn[3][0] }
	};
	gen_dn_dxy(revJ1, tmp, dn_dxy_1);
	//dn_dxy_2
	double tmp2[2][4] = {
		{derv_dN_dE[0][1],derv_dN_dE[1][1],derv_dN_dE[2][1],derv_dN_dE[3][1]},
		{derv_dN_dn[0][1], derv_dN_dn[1][1],derv_dN_dn[2][1],derv_dN_dn[3][1] }
	};
	gen_dn_dxy(revJ2, tmp2, dn_dxy_2);
	//dn_dxy_3
	double tmp3[2][4] = {
		{derv_dN_dE[0][2],derv_dN_dE[1][2],derv_dN_dE[2][2],derv_dN_dE[3][2]},
		{derv_dN_dn[0][2], derv_dN_dn[1][2],derv_dN_dn[2][2],derv_dN_dn[3][2] }
	};
	gen_dn_dxy(revJ3, tmp3, dn_dxy_3);
	//dn_dxy_4
	double tmp4[2][4] = {
		{derv_dN_dE[0][3],derv_dN_dE[1][3],derv_dN_dE[2][3],derv_dN_dE[3][3]},
		{derv_dN_dn[0][3], derv_dN_dn[1][3],derv_dN_dn[2][3],derv_dN_dn[3][3] }
	};
	gen_dn_dxy(revJ4, tmp4, dn_dxy_4);
	//dn_dxy_5
	double tmp5[2][4] = {
		{derv_dN_dE[0][4],derv_dN_dE[1][4],derv_dN_dE[2][4],derv_dN_dE[3][4]},
		{derv_dN_dn[0][4], derv_dN_dn[1][4],derv_dN_dn[2][4],derv_dN_dn[3][4] }
	};
	gen_dn_dxy(revJ5, tmp5, dn_dxy_5);
	//dn_dxy_6
	double tmp6[2][4] = {
		{derv_dN_dE[0][5],derv_dN_dE[1][5],derv_dN_dE[2][5],derv_dN_dE[3][5]},
		{derv_dN_dn[0][5],derv_dN_dn[1][5],derv_dN_dn[2][5],derv_dN_dn[3][5] }
	};
	gen_dn_dxy(revJ6, tmp6, dn_dxy_6);
	//dn_dxy_7
	double tmp7[2][4] = {
	{derv_dN_dE[0][6],derv_dN_dE[1][6],derv_dN_dE[2][6],derv_dN_dE[3][6]},
	{derv_dN_dn[0][6],derv_dN_dn[1][6],derv_dN_dn[2][6],derv_dN_dn[3][6] }
	};
	gen_dn_dxy(revJ7, tmp7, dn_dxy_7);
	//dn_dxy_8
	double tmp8[2][4] = {
	{derv_dN_dE[0][7],derv_dN_dE[1][7],derv_dN_dE[2][7],derv_dN_dE[3][7]},
	{derv_dN_dn[0][7],derv_dN_dn[1][7],derv_dN_dn[2][7],derv_dN_dn[3][7] }
	};
	gen_dn_dxy(revJ8, tmp8, dn_dxy_8);
	//dn_dxy_9
	double tmp9[2][4] = {
	{derv_dN_dE[0][8],derv_dN_dE[1][8],derv_dN_dE[2][8],derv_dN_dE[3][8]},
	{derv_dN_dn[0][8],derv_dN_dn[1][8],derv_dN_dn[2][8],derv_dN_dn[3][8] }
	};
	gen_dn_dxy(revJ9, tmp9, dn_dxy_9);
	//////////
	double d1_dx[4][4];
	double d1_dy[4][4];
	double tmpX1[4] = { dn_dxy_1[0][0],dn_dxy_1[0][1] ,dn_dxy_1[0][2] ,dn_dxy_1[0][3] };
	double tmpY1[4] = { dn_dxy_1[1][0],dn_dxy_1[1][1] ,dn_dxy_1[1][2] ,dn_dxy_1[1][3] };

	double d2_dx[4][4];
	double d2_dy[4][4];
	double tmpX2[4] = { dn_dxy_2[0][0],dn_dxy_2[0][1] ,dn_dxy_2[0][2] ,dn_dxy_2[0][3] };
	double tmpY2[4] = { dn_dxy_2[1][0],dn_dxy_2[1][1] ,dn_dxy_2[1][2] ,dn_dxy_2[1][3] };

	double d3_dx[4][4];
	double d3_dy[4][4];
	double tmpX3[4] = { dn_dxy_3[0][0],dn_dxy_3[0][1] ,dn_dxy_3[0][2] ,dn_dxy_3[0][3] };
	double tmpY3[4] = { dn_dxy_3[1][0],dn_dxy_3[1][1] ,dn_dxy_3[1][2] ,dn_dxy_3[1][3] };

	double d4_dx[4][4];
	double d4_dy[4][4];
	double tmpX4[4] = { dn_dxy_4[0][0],dn_dxy_4[0][1] ,dn_dxy_4[0][2] ,dn_dxy_4[0][3] };
	double tmpY4[4] = { dn_dxy_4[1][0],dn_dxy_4[1][1] ,dn_dxy_4[1][2] ,dn_dxy_4[1][3] };

	double d5_dx[4][4];
	double d5_dy[4][4];
	double tmpX5[4] = { dn_dxy_5[0][0],dn_dxy_5[0][1] ,dn_dxy_5[0][2] ,dn_dxy_5[0][3] };
	double tmpY5[4] = { dn_dxy_5[1][0],dn_dxy_5[1][1] ,dn_dxy_5[1][2] ,dn_dxy_5[1][3] };

	double d6_dx[4][4];
	double d6_dy[4][4];
	double tmpX6[4] = { dn_dxy_6[0][0],dn_dxy_6[0][1] ,dn_dxy_6[0][2] ,dn_dxy_6[0][3] };
	double tmpY6[4] = { dn_dxy_6[1][0],dn_dxy_6[1][1] ,dn_dxy_6[1][2] ,dn_dxy_6[1][3] };

	double d7_dx[4][4];
	double d7_dy[4][4];
	double tmpX7[4] = { dn_dxy_7[0][0],dn_dxy_7[0][1] ,dn_dxy_7[0][2] ,dn_dxy_7[0][3] };
	double tmpY7[4] = { dn_dxy_7[1][0],dn_dxy_7[1][1] ,dn_dxy_7[1][2] ,dn_dxy_7[1][3] };

	double d8_dx[4][4];
	double d8_dy[4][4];
	double tmpX8[4] = { dn_dxy_8[0][0],dn_dxy_8[0][1] ,dn_dxy_8[0][2] ,dn_dxy_8[0][3] };
	double tmpY8[4] = { dn_dxy_8[1][0],dn_dxy_8[1][1] ,dn_dxy_8[1][2] ,dn_dxy_8[1][3] };
	
	double d9_dx[4][4];
	double d9_dy[4][4];
	double tmpX9[4] = { dn_dxy_9[0][0],dn_dxy_9[0][1] ,dn_dxy_9[0][2] ,dn_dxy_9[0][3] };
	double tmpY9[4] = { dn_dxy_9[1][0],dn_dxy_9[1][1] ,dn_dxy_9[1][2] ,dn_dxy_9[1][3] };

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			d1_dx[j][i] = tmpX1[i] * tmpX1[j];
			d2_dx[j][i] = tmpX2[i] * tmpX2[j];
			d3_dx[j][i] = tmpX3[i] * tmpX3[j];
			d4_dx[j][i] = tmpX4[i] * tmpX4[j];
			d5_dx[j][i] = tmpX5[i] * tmpX5[j];
			d6_dx[j][i] = tmpX6[i] * tmpX6[j];
			d7_dx[j][i] = tmpX7[i] * tmpX7[j];
			d8_dx[j][i] = tmpX8[i] * tmpX8[j];
			d9_dx[j][i] = tmpX9[i] * tmpX9[j];

			d1_dy[j][i] = tmpY1[i] * tmpY1[j];
			d2_dy[j][i] = tmpY2[i] * tmpY2[j];
			d3_dy[j][i] = tmpY3[i] * tmpY3[j];
			d4_dy[j][i] = tmpY4[i] * tmpY4[j];
			d5_dy[j][i] = tmpY5[i] * tmpY5[j];
			d6_dy[j][i] = tmpY6[i] * tmpY6[j];
			d7_dy[j][i] = tmpY7[i] * tmpY7[j];
			d8_dy[j][i] = tmpY8[i] * tmpY8[j];
			d9_dy[j][i] = tmpY9[i] * tmpY9[j];
		}
	}
	/*
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << d1_dx[j][i] << " ";
		}
		std::cout << std::endl;
	}
	*/

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			d1_dx[j][i] += d1_dy[j][i];
			d2_dx[j][i] += d2_dy[j][i];
			d3_dx[j][i] += d3_dy[j][i];
			d4_dx[j][i] += d4_dy[j][i];
			d5_dx[j][i] += d5_dy[j][i];
			d6_dx[j][i] += d6_dy[j][i];
			d7_dx[j][i] += d7_dy[j][i];
			d8_dx[j][i] += d8_dy[j][i];
			d9_dx[j][i] += d9_dy[j][i];
		}
	}
	//mnozenie razy k
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			d1_dx[j][i] *= this->k_t;
			d2_dx[j][i] *= this->k_t;
			d3_dx[j][i] *= this->k_t;
			d4_dx[j][i] *= this->k_t;
			d5_dx[j][i] *= this->k_t;
			d6_dx[j][i] *= this->k_t;
			d7_dx[j][i] *= this->k_t;
			d8_dx[j][i] *= this->k_t;
			d9_dx[j][i] *= this->k_t;
		}
	}
	//* det
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			d1_dx[j][i] *= detJ1 * w[0] * w[0];
			d2_dx[j][i] *= detJ2 * w[1] * w[0];
			d3_dx[j][i] *= detJ3 * w[2] * w[0];
			d4_dx[j][i] *= detJ4 * w[0] * w[1];
			d5_dx[j][i] *= detJ5 * w[1] * w[1];
			d6_dx[j][i] *= detJ6 * w[2] * w[1];
			d7_dx[j][i] *= detJ7 * w[0] * w[2];
			d8_dx[j][i] *= detJ8 * w[1] * w[2];
			d9_dx[j][i] *= detJ9 * w[2] * w[2];
		}
	}
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			H[j][i] += d1_dx[j][i] ; 
			H[j][i] += d2_dx[j][i] ;
			H[j][i] += d3_dx[j][i] ;
			H[j][i] += d4_dx[j][i] ;
			H[j][i] += d5_dx[j][i] ;
			H[j][i] += d6_dx[j][i] ;
			H[j][i] += d7_dx[j][i] ;
			H[j][i] += d8_dx[j][i] ;
			H[j][i] += d9_dx[j][i] ;
		}
	}

}




#endif // !ELEM_UNIWERSAL_3_Point
