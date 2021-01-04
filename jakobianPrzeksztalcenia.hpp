#ifndef JAKOBIAN_PRZEKSZTALCENIA
#define JAKOBIAN_PRZEKSZTALCENIA
#include <iostream>

class ElemUniwersal4_2point {
private:
	//								1 / sqrt(3) ( punkty calkowania )
	double E[4] = { -0.5773502692, 0.5773502692, 0.5773502692, -0.5773502692 }; //like X
	double n[4] = { -0.5773502692, -0.5773502692, 0.5773502692, 0.5773502692 }; //like Y
	//
	//jakobiany dla kazdego z 4 ech punktow
	double J1[2][2];
	double J2[2][2];
	double J3[2][2];
	double J4[2][2];

	double revJ1[2][2];
	double revJ2[2][2];
	double revJ3[2][2];
	double revJ4[2][2];

	double detJ1;
	double detJ2;
	double detJ3;
	double detJ4;

	//dN/dE
	double derv_dN_dE[4][4];
	//dN/dn
	double derv_dN_dn[4][4];

	double dn_dxy_1[2][4];
	double dn_dxy_2[2][4];
	double dn_dxy_3[2][4];
	double dn_dxy_4[2][4];

	double H[4][4];
	double C[4][4];
	double H_BC[4][4];
	double N[4][4];

	void cTmpTimesDet(double** cTmp, int detNumber);
	void addToC(double** cTmp);
	void H_BCTimesDet(int wall, double det);
	void H_BCTimesAlfa();

	double P[4];
	void PTimesDet(int wall, double det);
	void PTimesAlfaT_inf();

	double k_t = 25;
	double alfa = 300;
	double t_inf = 1200;
	double deltaTau = 50;

	void genJ(double* xArr, double* yArr);
	void genRevJ();
	void gen_dn_dxy(double revJ[2][2], double tmp[2][4], double dn_dxy[2][4]);
public:
	void genJacobians(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);
	void jClear(int numOfJ);
	void jPrint(int numOfJ);
	void derv_dN_dEPrint();
	void derv_dN_dnPrint();
	void revJPrint(int numOfJ);
	void dn_xyPrint(int numOfdx_dx);
	void hPrint();
	void genC(double c, double ro);
	void CClear();
	void CPrint();
	void genH_BC(double x1, double y1, int BC1, double x2, double y2, int BC2, double x3, double y3, int BC3, double x4, double y4, int BC4);
	void H_BCClear();
	void H_BCPrint();
	void genP(double x1, double y1, int BC1, double x2, double y2, int BC2, double x3, double y3, int BC3, double x4, double y4, int BC4);
	void PClear();
	void PPrint();

	double getT_inf() { return t_inf; }
	double getAlfa() { return alfa; }
	double getDeltaTau() { return deltaTau; }
	
	void setT_inf(double T_inf) { t_inf = T_inf; }
	void setAlfa(double alfa) { this->alfa = alfa; }
	void setDeltaTau(double deltaTau) { this->deltaTau = deltaTau; }

	void addH_BC2H();

	double** getH();
	double** getC();
	double* getP();
	double** getH_BC();

	ElemUniwersal4_2point() = default;
};

void ElemUniwersal4_2point::addH_BC2H() {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			H[i][j] += H_BC[i][j];
		}
	}
}

void ElemUniwersal4_2point::PClear() {
	for (int i = 0; i < 4; ++i) {
		P[i] = 0;
	}
}

void ElemUniwersal4_2point::PPrint() {
	std::cout << std::endl << "P: " << std::endl;
	for (int i = 0; i < 4; ++i) {
		std::cout << P[i] << "		";
	}
}

void ElemUniwersal4_2point::PTimesDet(int wall, double det) {
	if (wall == 0) {
		P[0] *= det;
		P[1] *= det;
	}
	if (wall == 1) {
		P[1] *= det;
		P[2] *= det;
	}
	if (wall == 2) {
		P[2] *= det;
		P[3] *= det;
	}
	if (wall == 3) {
		P[3] *= det;
		P[0] *= det;
	}
}

void ElemUniwersal4_2point::PTimesAlfaT_inf() {
	for (int i = 0; i < 4; ++i) {
		P[i] *= -alfa;
		P[i] *= t_inf;
	}
}

void ElemUniwersal4_2point::genP(double x1, double y1, int BC1, double x2, double y2, int BC2, double x3, double y3, int BC3, double x4, double y4, int BC4) {
	PClear();
	double N[2][2];
	double L;
	for (int i = 0; i < 2; ++i) {
		N[i][0] = (1 - E[i]) / 2;
		N[i][1] = (1 + E[i]) / 2;
	}

	//dolny bok
	if (BC1 == 1 && BC2 == 1) {
		L = (x2 - x1) / 2;
		for (int i = 0; i < 2; ++i) {
			P[0] += N[i][0] * L;
			P[1] += N[i][1] * L;
		}
		//PTimesDet(0, L);
	}
	//prawy bok
	if (BC2 == 1 && BC3 == 1) {
		L = (y3 - y2) / 2;
		for (int i = 0; i < 2; ++i) {
			P[1] += N[i][0] * L;
			P[2] += N[i][1] * L;
		}
		//PTimesDet(1, L);
	}
	//gorny bok
	if (BC3 == 1 && BC4 == 1) {
		L = (x3 - x4) / 2;
		for (int i = 0; i < 2; ++i) {
			P[2] += N[i][0] * L;
			P[3] += N[i][1] * L;
		}
		//PTimesDet(2, L);
	}
	//lewy bok
	if (BC4 == 1 && BC1 == 1) {
		L = (y4 - y1) / 2;
		for (int i = 0; i < 2; ++i) {
			P[3] += N[i][0] * L;
			P[0] += N[i][1] * L;
		}
		//PTimesDet(3, L);
	}

	PTimesAlfaT_inf();
}

void ElemUniwersal4_2point::H_BCTimesAlfa() {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			H_BC[i][j] *= 25;//alfa;
		}
	}
}

void ElemUniwersal4_2point::H_BCTimesDet(int wall, double det) {
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

void ElemUniwersal4_2point::H_BCPrint() {
	std::cout << std::endl << "H_BC" << std::endl;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << H_BC[i][j] << "		";
		}
		std::cout << std::endl;
	}
}

void ElemUniwersal4_2point::H_BCClear() {
	
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			H_BC[i][j] = 0;
		}
	}
	
}

void ElemUniwersal4_2point::genH_BC(double x1, double y1, int BC1, double x2, double y2, int BC2, double x3, double y3, int BC3, double x4, double y4, int BC4) {
	H_BCClear();
	double L;
	double N1, N2;
	N1 = (1 - E[0]) / 2;
	N2 = (1 + E[0]) / 2;
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
		tmp1[0][0] += N1 * N1;
		tmp1[0][1] += N1 * N2;
		tmp1[1][0] += N1 * N2;
		tmp1[1][1] += N2 * N2;

		N1 = (1 - E[1]) / 2;
		N2 = (1 + E[1]) / 2;
		tmp1[0][0] += N1 * N1;
		tmp1[0][1] += N1 * N2;
		tmp1[1][0] += N1 * N2;
		tmp1[1][1] += N2 * N2;
		
		tmp1[0][0] *= L;
		tmp1[0][1] *= L;
		tmp1[1][0] *= L;
		tmp1[1][1] *= L;
		//H_BCTimesDet(0, L);
	}
	//prawy bok
	if (BC2 == 1 && BC3 == 1) {
		N1 = (1 - E[0]) / 2;
		N2 = (1 + E[0]) / 2;
		L = (y3 - y2) / 2;
		tmp2[0][0] += N1 * N1;
		tmp2[0][1] += N1 * N2;
		tmp2[1][0] += N1 * N2;
		tmp2[1][1] += N2 * N2;

		N1 = (1 - E[1]) / 2;
		N2 = (1 + E[1]) / 2;
		tmp2[0][0] += N1 * N1;
		tmp2[0][1] += N1 * N2;
		tmp2[1][0] += N1 * N2;
		tmp2[1][1] += N2 * N2;
		
		tmp2[0][0] *= L;
		tmp2[0][1] *= L; 
		tmp2[1][0] *= L;
		tmp2[1][1] *= L;
		//H_BCTimesDet(1, L);
	}
	//gorny bok
	if (BC3 == 1 && BC4 == 1) {
		N1 = (1 - E[0]) / 2;
		N2 = (1 + E[0]) / 2;
		L = (x3 - x4) / 2;
		tmp3[0][0] += N1 * N1;
		tmp3[0][1] += N1 * N2;
		tmp3[1][0] += N1 * N2;
		tmp3[1][1] += N2 * N2;
	
		N1 = (1 - E[1]) / 2;
		N2 = (1 + E[1]) / 2;
		tmp3[0][0] += N1 * N1;
		tmp3[0][1] += N1 * N2;
		tmp3[1][0] += N1 * N2;
		tmp3[1][1] += N2 * N2;
		
		tmp3[0][0] *= L;
		tmp3[0][1] *= L;
		tmp3[1][0] *= L;
		tmp3[1][1] *= L;
		//H_BCTimesDet(2, L);
	}
	//lewy bok
	if (BC4 == 1 && BC1 == 1) {
		N1 = (1 - E[0]) / 2;
		N2 = (1 + E[0]) / 2;
		L = (y4 - y1) / 2;
		tmp4[0][0] += N1 * N1;
		tmp4[0][1] += N1 * N2;
		tmp4[1][0] += N1 * N2;
		tmp4[1][1] += N2 * N2;

		N1 = (1 - E[1]) / 2;
		N2 = (1 + E[1]) / 2;
		tmp4[0][0] += N1 * N1;
		tmp4[0][1] += N1 * N2;
		tmp4[1][0] += N1 * N2;
		tmp4[1][1] += N2 * N2;
		
		tmp4[0][0] *= L;
		tmp4[0][1] *= L;
		tmp4[1][0] *= L;
		tmp4[1][1] *= L;
		//H_BCTimesDet(3, L);
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

void ElemUniwersal4_2point::CPrint() {
	std::cout << std::endl << "C:" << std::endl;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << C[i][j] << "			";
		}
		std::cout << std::endl;
	}
}

void ElemUniwersal4_2point::CClear() {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			C[i][j] = 0;
		}
	}
}

void ElemUniwersal4_2point::addToC(double** cTmp) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			C[i][j] += cTmp[i][j];
		}
	}
}


double* ElemUniwersal4_2point::getP() {
	double* tmp = new double[4];
	for (int i = 0; i < 4; ++i) {
		tmp[i] = P[i];
	}
	return tmp;
}

void ElemUniwersal4_2point::cTmpTimesDet(double** cTmp, int detNumber) {
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			switch (detNumber) {
			case 1:
				cTmp[i][j] *= detJ1;
				break;
			case 2:
				cTmp[i][j] *= detJ2;
				break;
			case 3:
				cTmp[i][j] *= detJ3;
				break;
			case 4:
				cTmp[i][j] *= detJ4;
				break;
			default:
				break;
			}
		}
	}
}

void ElemUniwersal4_2point::genC(double c, double ro) {

	CClear();

	for (int i = 0; i < 4; ++i) {
		N[0][i] = 0.25 * (1 - E[i]) * (1 - n[i]);
		N[1][i] = 0.25 * (1 + E[i]) * (1 - n[i]);
		N[2][i] = 0.25 * (1 + E[i]) * (1 + n[i]);
		N[3][i] = 0.25 * (1 - E[i]) * (1 + n[i]);
	}
	double** cTmp = new double*[4];
	for (int i = 0; i < 4; ++i) {
		cTmp[i] = new double[4];
	}


	//for each row
	for (int i = 0; i < 4; ++i) {
		//clear cTmp
		for (int k = 0; k < 4; ++k) {
			for (int j = 0; j < 4; ++j) {
				cTmp[j][k] = 0;
			}
		}

		for (int a = 0; a < 4; ++a) {
			for (int b = 0; b < 4; ++b) {
				cTmp[b][a] = N[a][i] * N[b][i];// * c * ro;
			}
		}
		/*for (int a = 0; a < 4; ++a) {
			cTmp[a][0] = N[0][i] * N[0][a];
		}*/
		//cTmpTimesDet(cTmp, i);
		
		addToC(cTmp);

		/*std::cout << i << "wart c" << std::endl;
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				std::cout << cTmp[i][j] << "		";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;*/
	}
	
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			C[i][j] *= c * ro * detJ1;
		}
	}
	
}
double** ElemUniwersal4_2point::getC() {
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

double** ElemUniwersal4_2point::getH_BC() {
	double **tmp = new double*[4];
	for (int i = 0; i < 4; ++i) {
		tmp[i] = new double[4];
	}
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			tmp[j][i] = H_BC[j][i];
		}
	}
	return tmp;
}

double** ElemUniwersal4_2point::getH() {
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

void ElemUniwersal4_2point::hPrint( ) {
	std::cout << std::endl << "H" << std::endl;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << H[j][i] << " ";
		}
		std::cout << std::endl;
	}
	
}


void ElemUniwersal4_2point::revJPrint(int numOfJ) {
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
		break;

	}
}

void ElemUniwersal4_2point::dn_xyPrint(int numOfdx_dy) {
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

	}

}

void ElemUniwersal4_2point::jPrint(int numOfJ) {
	std::cout << std::endl << "J" << numOfJ << std::endl;
	switch (numOfJ)
	{
	case 1:
		for (int i = 0; i < 2; ++i) {
			std::cout << "| " << J1[i][0]  << "      ";
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
	default:
		break;
	}
}

void ElemUniwersal4_2point::jClear(int numOfJ) {
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
	default:
		break;
	}
}

void ElemUniwersal4_2point::derv_dN_dEPrint() {
	std::cout << std::endl;
	std::cout << "derv_dN/dE" << std::endl;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << derv_dN_dE[j][i] << " ";
		}
		std::cout << std::endl;
	}


}

void ElemUniwersal4_2point::derv_dN_dnPrint() {
	std::cout << std::endl;
	std::cout << "derv_dN/dn" << std::endl;
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			std::cout << derv_dN_dn[j][i] << " ";
		}
		std::cout << std::endl;
	}
}

void ElemUniwersal4_2point::genJ(double* xArr, double* yArr) {
	for (int i = 0; i < 4; ++i) {
		//	dX/dE
		J1[0][0] += derv_dN_dE[i][0] * xArr[i];
		//	dY/dE
		J1[1][0] += derv_dN_dE[i][0] * yArr[i];
		//	dX/dn
		J1[0][1] += derv_dN_dn[i][0] * xArr[i];
		//	dY/dn
		J1[1][1] += derv_dN_dn[i][0] * yArr[i];
	
	/////
	//J2
	
		//	dX/dE
		J2[0][0] += derv_dN_dE[i][1] * xArr[i];
		//	dY/dE
		J2[1][0] += derv_dN_dE[i][1] * yArr[i];
		//	dX/dn
		J2[0][1] += derv_dN_dn[i][1] * xArr[i];
		//	dY/dn
		J2[1][1] += derv_dN_dn[i][1] * yArr[i];
	
	/////
	//J3

		//	dX/dE
		J3[0][0] += derv_dN_dE[i][2] * xArr[i];
		//	dY/dE
		J3[1][0] += derv_dN_dE[i][2] * yArr[i];
		//	dX/dn
		J3[0][1] += derv_dN_dn[i][2] * xArr[i];
		//	dY/dn
		J3[1][1] += derv_dN_dn[i][2] * yArr[i];
	
	/////
	//J2

		//	dX/dE
		J4[0][0] += derv_dN_dE[i][3] * xArr[i];
		//	dY/dE
		J4[1][0] += derv_dN_dE[i][3] * yArr[i];
		//	dX/dn
		J4[0][1] += derv_dN_dn[i][3] * xArr[i];
		//	dY/dn
		J4[1][1] += derv_dN_dn[i][3] * yArr[i];
	}
}

void ElemUniwersal4_2point::genRevJ() {
		detJ1 = (J1[0][0] * J1[1][1]) - (J1[0][1] * J1[1][0]);
		revJ1[0][0] = J1[1][1] * (1 / detJ1);
		revJ1[1][0] = -J1[1][0] * (1 / detJ1);
		revJ1[0][1] = -J1[0][1] * (1 / detJ1);
		revJ1[1][1] = J1[0][0] * (1 / detJ1);
	
	
		detJ2 = (J2[0][0] * J2[1][1]) - (J2[0][1] * J2[1][0]);
		revJ2[0][0] = J2[1][1] * (1 / detJ2);
		revJ2[1][0] = -J2[1][0] * (1 / detJ2);
		revJ2[0][1] = -J2[0][1] * (1 / detJ2);
		revJ2[1][1] = J2[0][0] * (1 / detJ2);
	
	
		detJ3 = (J3[0][0] * J3[1][1]) - (J3[0][1] * J3[1][0]);
		revJ3[0][0] = J3[1][1] * (1 / detJ3);
		revJ3[1][0] = -J3[1][0] * (1 / detJ3);
		revJ3[0][1] = -J3[0][1] * (1 / detJ3);
		revJ3[1][1] = J3[0][0] * (1 / detJ3);
	
	
		detJ4 = (J4[0][0] * J4[1][1]) - (J4[0][1] * J4[1][0]);
		revJ4[0][0] = J4[1][1] * (1 / detJ4);
		revJ4[1][0] = -J4[1][0] * (1 / detJ4);
		revJ4[0][1] = -J4[0][1] * (1 / detJ4);
		revJ4[1][1] = J4[0][0] * (1 / detJ4);			
}

void ElemUniwersal4_2point::gen_dn_dxy(double revJ[2][2], double tmp[2][4], double dn_dxy[2][4]) {
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

void ElemUniwersal4_2point::genJacobians(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4) {
	double xArr[4] = { x1,x2,x3,x4 };
	double yArr[4] = { y1,y2,y3,y4 };

	for (int i = 0; i < 4; ++i) {
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
	//J
	jClear(1);
	jClear(2);
	jClear(3);
	jClear(4);
	genJ(xArr, yArr);
	
	//revJ1
	genRevJ();

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

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			d1_dx[j][i] = tmpX1[i] * tmpX1[j];
			d2_dx[j][i] = tmpX2[i] * tmpX2[j];
			d3_dx[j][i] = tmpX3[i] * tmpX3[j];
			d4_dx[j][i] = tmpX4[i] * tmpX4[j];

			d1_dy[j][i] = tmpY1[i] * tmpY1[j];
			d2_dy[j][i] = tmpY2[i] * tmpY2[j];
			d3_dy[j][i] = tmpY3[i] * tmpY3[j];
			d4_dy[j][i] = tmpY4[i] * tmpY4[j];
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
		}
	}
	//mnozenie razy k i det
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			d1_dx[j][i] *= this->k_t * detJ1;
			d2_dx[j][i] *= this->k_t * detJ2;
			d3_dx[j][i] *= this->k_t * detJ3;
			d4_dx[j][i] *= this->k_t * detJ4;
		}
	}

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			H[j][i] += d1_dx[j][i]; 
			H[j][i] += d2_dx[j][i];
			H[j][i] += d3_dx[j][i];
			H[j][i] += d4_dx[j][i];
		}
	}

}

#endif // !JAKOBIAN_PRZEKSZTALCENIA
