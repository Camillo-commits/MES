#include <iostream>
#include<conio.h>
#include <string>
#include<fstream>

#include"netGen.hpp"
#include"calkNum.hpp"
#include"jakobianPrzeksztalcenia.hpp"
#include"ElementUniwersal_3Point.h"
#include"Jacobi.h"

double fun(double x, double y) {
	return -2*x*x*y + 2*x*y + 4;
}
double fun2(double x, double y) {
	return -5 * x*x*y + 2*x*y*y + 10;
}

int main() {

	GlobalData* gD = new GlobalData();
	gD->setAlfa(300);
	gD->setT_inf(1200);
	gD->setDeltaTau(50);
	
	Node* n = new Node [gD->getNN()];
	Element* e = new Element [gD->getNE()];

	for (int i = 0; i < gD->getNN();++i) {
		n[i] = Node();
	}
	for (int i = 0; i < gD->getNE(); ++i) {
		e[i] = Element();
	}

	gD->genArrays(e, n,100);
	for (int i = 0; i < gD->getNH(); ++i) {
		for (int j = 0; j < gD->getNW(); ++j) {
			std::cout << n[(i + (gD->getNW() * j))].getX() << "," << n[i].getY() << " ";
		}
		std::cout << std::endl;
	}

	std::cout << "nodes: " << std::endl;
	for (int i = 0; i < gD->getNN(); ++i) {
		std::cout << n[i].getX() << " , " << n[i].getY() << std::endl;
	}
	

	/*
	for (int elem = 0; elem < gD->getNE() ; ++elem) {
		std::cout << "elem " << elem + 1<< std::endl;
		for (int node = 0; node < 4; ++node) {
			std::cout << "id = " << e[elem].getId(node)->getGlobalId() << "		x = " << e[elem].getId(node)->getX() << "		y = " << e[elem].getId(node)->getY() << std::endl;
		}
		std::cout << std::endl;
	}
	*/
	/*
	//test 3 punktowy 2 punktowy
	ElemUniwersal4_2point* elemUniwersal = new ElemUniwersal4_2point();
	elemUniwersal->genJacobians(0,0,4,0,4,6,0,6);
	elemUniwersal->hPrint();
	*/
	/*
	ElemUniwersal4_2point* elemUniwersal2Point = new ElemUniwersal4_2point();
	//elemUniwersal3Point->genJacobians(0, 0, 4, 0, 4, 6, 0, 6);
	elemUniwersal2Point->genJacobians(0, 0, 0.2, 0, 0.2, 0.2, 0, 0.2);
	//elemUniwersal3Point->hPrint();

	elemUniwersal2Point->genC(700, 7800);
	elemUniwersal2Point->CPrint();

	ElemUniwersal4_3point* elemUniwersal3Point = new ElemUniwersal4_3point();
	//elemUniwersal3Point->genJacobians(0, 0, 4, 0, 4, 6, 0, 6);
	elemUniwersal3Point->genJacobians(0, 0, 0.2, 0, 0.2, 0.2, 0, 0.2);
	//elemUniwersal3Point->hPrint();

	elemUniwersal3Point->genC(700, 7800);
	elemUniwersal3Point->CPrint();
	*/
	
	//////////////////////////////

	//////////////////////////////
	/*
	ElemUniwersal4_2point* elemUniwersal2 = new ElemUniwersal4_2point();
	elemUniwersal2->genJacobians(1, 1, 5, 1, 5, 7, 1, 7);
	elemUniwersal2->hPrint();
	*/
	//////////////////////////////
	
	////ustawienie zmiennej BC na œcianach

	

	/////////////////////////////////////

	
	SOE* soe = new SOE(gD->getNN());
	//sets starting temp info
	soe->setT(n, gD->getNN());
	soe->sym(e, gD->getNE(), gD->getNN(), 3, gD->getAlfa(), gD->getT_inf(), gD->getDeltaTau(), 500);
	
	/*soe->genHCP_global(e, gD->getNE(),gD->getNN(),2, gD->getAlfa(), gD->getT_inf(),gD->getDeltaTau());
	soe->H_Print(gD->getNN());
	soe->C_Print(gD->getNN());

	soe->genHCP_global(e, gD->getNE(), gD->getNN(), 3, gD->getAlfa(), gD->getT_inf(),gD->getDeltaTau());
	soe->H_Print(gD->getNN());
	soe->C_Print(gD->getNN());
	//test macierzy hbc
	/*
	ElemUniwersal4_2point* elemUniwersal2 = new ElemUniwersal4_2point();
	ElemUniwersal4_3point* elemUniwersal3 = new ElemUniwersal4_3point();

	//lewa
	elemUniwersal2->genH_BC(0, 0, 1, 0.03333, 0, 0, 0.03333, 0.03333, 0, 0, 0.033333, 1);
	elemUniwersal2->H_BCPrint();
	elemUniwersal3->genH_BC(0, 0, 1, 0.03333, 0, 0, 0.03333, 0.03333, 0, 0, 0.033333, 1);
	elemUniwersal3->H_BCPrint();

	/* //dol
	elemUniwersal2->genH_BC(0, 0, 1, 0.03333, 0, 1, 0.03333, 0.03333, 0, 0, 0.033333, 0);
	elemUniwersal2->H_BCPrint();
	elemUniwersal3->genH_BC(0, 0, 1, 0.03333, 0, 1, 0.03333, 0.03333, 0, 0, 0.033333, 0);
	elemUniwersal3->H_BCPrint();

	//prawa
	elemUniwersal2->genH_BC(0, 0, 0, 0.03333, 0, 1, 0.03333, 0.03333, 1, 0, 0.033333, 0);
	elemUniwersal2->H_BCPrint();
	elemUniwersal3->genH_BC(0, 0, 0, 0.03333, 0, 1, 0.03333, 0.03333, 1, 0, 0.033333, 0);
	elemUniwersal3->H_BCPrint();
	*/
	//gora
	/*
	elemUniwersal2->genH_BC(0, 0, 0, 0.03333, 0, 0, 0.03333, 0.03333, 1, 0, 0.033333, 1);
	elemUniwersal2->H_BCPrint();

	elemUniwersal3->genH_BC(0, 0, 0, 0.03333, 0, 0, 0.03333, 0.03333, 1, 0, 0.033333, 1);
	elemUniwersal3->H_BCPrint();
	////////////////
	*/

	/*
	////test calki
	std::cout << "calka" << std::endl;
	std::cout << calka2DStandarized2Points(fun) <<std::endl;
	std::cout << calka2DStandarized3Points(fun2) << std::endl;
	////////////
	*/
	/*
	ElemUniwersal4* elemUniwersal4 = new ElemUniwersal4();
	//elemUniwersal4 -> genJacobians(0, 0, 4, 0, 4, 4, 0, 4);
	elemUniwersal4->genJacobians(0, 0, 4, 0, 4, 6, 0, 6);
	elemUniwersal4 -> jPrint(1);
	elemUniwersal4->jPrint(2);
	elemUniwersal4->jPrint(3);
	elemUniwersal4->jPrint(4);

	elemUniwersal4->derv_dN_dEPrint();
	elemUniwersal4->derv_dN_dnPrint();

	elemUniwersal4->revJPrint(1);
	elemUniwersal4->revJPrint(2);
	elemUniwersal4->revJPrint(3);
	elemUniwersal4->revJPrint(4);

	elemUniwersal4->dn_xyPrint(1);
	elemUniwersal4->dn_xyPrint(2);
	elemUniwersal4->dn_xyPrint(3);
	elemUniwersal4->dn_xyPrint(4);

	elemUniwersal4->hPrint();
	*/
	
	
	/*double** X = new double*[2];
	for (int i = 0; i < 2; ++i) {
		X[i] = new double[2];
	}
	X[0][0] = 1;
	X[0][1] = 1;
	X[1][0] = 4;
	X[1][1] = 0;
	double* B = new double[2];
	B[0] = 2;
	B[1] = 4;
	
	B = resolveSOEJacobi(X, B, 2);
	std::cout << B[0] << "				" << B[1];
	*/
	_getch();
	return 0;
}