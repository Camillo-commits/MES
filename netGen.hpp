#ifndef NET_GEN
#define NET_GEN


//#include"netGen.hpp"
#include<fstream>
#include<string>
#include<cmath>
#include <iostream>
#include"jakobianPrzeksztalcenia.hpp"
#include "ElementUniwersal_3Point.h"
#include"Jacobi.h"


class Node {
private:
	double x, y;
	int globalId;
	int BC;
	double t;
public:
	Node(double X, double Y, int GLOBALID) : x(X), y(Y), globalId(GLOBALID) { BC = 0; };
	Node(double X, double Y, int GLOBALID, int BC) : x(X), y(Y), globalId(GLOBALID), BC(BC) {};
	Node() :x(0), y(0), globalId(-1) { BC = 0; };

	double getX() { return x; }
	double getY() { return y; }
	int getGlobalId() { return globalId; }
	int getBC() { return BC; }
	double getT() { return t; }

	void setGlobalId(int id) { globalId = id; }
	void setX(double X) { x = X; }
	void setY(double Y) { y = Y; }
	void setBC(int BC) { this->BC = BC; }
	void setT(int T) { t = T; }
};

class Element {
private:
	Node* id[4];
	double H_lok[4][4];
	double C_lok[4][4];
	double P_lok[4];
	//ElemUniwersal4* elemUniwersal4;
public:
	Element(Node* ID0, Node* ID1, Node* ID2, Node* ID3);
	Element() { id[0] = new Node(); id[1] = new Node(); id[2] = new Node(); id[3] = new Node(); };

	Node* getId(int i) { return id[i]; }
	void setId(Node ID, int i) {
		id[i]->setX(ID.getX());
		id[i]->setY(ID.getY());
		id[i]->setGlobalId(ID.getGlobalId());
	}
	double getH_lok(int wiersz, int kolumna) {
		return H_lok[wiersz][kolumna];
	}
	void setH_lok(double **H, int numOfWiersz, int numOfKolumna);

	double getC_lok(int wiersz, int kolumna) {
		return C_lok[wiersz][kolumna];
	}
	void setC_lok(double **C, int numOfWiersz, int numOfKolumna);
	double getP_lok(int kolumna) {
		return P_lok[kolumna];
	}
	void setP_lok(double *P, int numOfColumns);
	//ElemUniwersal4* getElemUniwersal4() { return elemUniwersal4; }
	//void setElemUniwersal4(ElemUniwersal4* elem) { elemUniwersal4 = elem; }
};

void Element::setP_lok(double *P, int numOfColumns) {
	for (int kolumna = 0; kolumna < numOfColumns; ++kolumna) {
		this->P_lok[kolumna] = P[kolumna];
	}
}

void Element::setC_lok(double **C, int numOfWiersz, int numOfKolumna) {
	for (int kolumna = 0; kolumna < numOfKolumna; ++kolumna) {
		for (int wiersz = 0; wiersz < numOfWiersz; ++wiersz) {
			this->C_lok[wiersz][kolumna] = C[wiersz][kolumna];
		}
	}
}

void Element::setH_lok(double **H, int numOfWiersz, int numOfKolumna) {
	for (int kolumna = 0; kolumna < numOfKolumna; ++kolumna) {
		for (int wiersz = 0; wiersz < numOfWiersz; ++wiersz) {
			this->H_lok[wiersz][kolumna] = H[wiersz][kolumna];
		}
	}

}

Element::Element(Node* ID0, Node* ID1, Node* ID2, Node* ID3) {
	id[0] = ID0; id[1] = ID1; id[2] = ID2; id[3] = ID3; 
	//elemUniwersal4 = new ElemUniwersal4();
	//elemUniwersal4->genJacobians(ID0->getX(), ID0->getY(), ID1->getX(), ID1->getY(), ID2->getX(), ID2->getY(), ID3->getX(), ID3->getY());
}

class GlobalData {
private:
	//H -> Height, W -> Width, nH -> number of elements from top to bottom, nW -> number of elements from left to right, nE -> number of elements, nN -> number of nodes
	double H, W;
	unsigned int nH, nW;
	unsigned int nE, nN;
	double t_inf;
	double alfa;
	double deltaTau;
	double k_t;
	void setBC(Element* E);
public:
	GlobalData();

	double getH() { return H; }
	double getW() { return W; }
	unsigned int getNH() { return nH; }
	unsigned int getNW() { return nW; }
	unsigned int getNE() { return nE; }
	unsigned int getNN() { return nN; }
	double getAlfa() { return alfa; }
	double getT_inf() { return t_inf; }
	double getDeltaTau() { return deltaTau; }
	double getK_t() { return k_t; }

	void setH(double h) { H = h; }
	void setW(double w) { W = w; }
	void setNH(unsigned int nh) { nH = nh; }
	void setNW(unsigned int nw) { nW = nw; }
	void setNE(unsigned int ne) { nE = ne; }
	void setNN(unsigned int nn) { nN = nn; }
	void setAlfa(double alfa) { this->alfa = alfa; }
	void setT_inf(double t_inf) { this->t_inf = t_inf; }
	void setDeltaTau(double deltaTau) { this->deltaTau = deltaTau; }
	void setK_t(double k_t) { this->k_t = k_t; }

	void genArrays(Element* E, Node* N, double t0);
};

class SOE {
private:
	double** H_global;
	double** C_global;
	double* P_global;
	double* t;		//temperatura

	void modify(Element* E, double deltaTau, int numOfElem, int numOfNodes);
public:
	SOE();
	SOE(int sizeOfH);
	double getElemFromH(int row, int column);
	void setElemFromH(int row, int column, double value);
	void addElemToH(int row, int column, double value);
	void genHCP_global(Element* e, int numOfElem,int numOfNodes,int numOfIntegralPoints, double alfa, double t_inf,double deltaTau);
	void H_Clear(int numOfElem);
	void H_Print(int sizeOfH);
	void C_Clear(int numOfNodes);
	void C_Print(int sizeOfC);
	void P_Clear(int numOfNodes);
	void P_Print(int sizeOfP);
	double* getT(int numOfNodes);
	void setT(Node* n,int numOfNodes);
	void sym(Element* e, int numOfElem, int numOfNodes, int numOfIntegralPoints, double alfa, double t_inf, double deltaTau,double maxT);
	void t_Print(int sizeOfT);
};

void SOE::sym(Element* e, int numOfElem, int numOfNodes, int numOfIntegralPoints, double alfa, double t_inf, double deltaTau,double maxT) {
	std::cout << "t0: "; t_Print(numOfNodes);
	for (int i = 0; i < maxT / deltaTau; ++i) {
		genHCP_global(e, numOfElem, numOfNodes, numOfIntegralPoints, alfa, t_inf,deltaTau);
		t = resolveSOEJacobi(H_global, P_global, numOfNodes);
		std::cout << "t" << i;
		t_Print(numOfNodes);
	}
}

void SOE::setT(Node* n, int numOfNodes) {
	for (int i = 0; i < numOfNodes; ++i) {
		t[i] = n[i].getT();
	}
}

double* SOE::getT(int numOfNodes) {
	double* tmp = new double[numOfNodes];
	for (int i = 0; i < numOfNodes; ++i) {
		tmp[i] = t[i];
	}
	return tmp;
}

void SOE::C_Print(int sizeOfC) {
	std::cout << "C global:" << std::endl;
	for (int i = 0; i < sizeOfC; ++i) {
		for (int j = 0; j < sizeOfC; ++j) {
			std::cout << C_global[j][i] << "	";
		}
		std::cout << std::endl;
	}
}

void SOE::C_Clear(int numOfNodes) {
	for (int i = 0; i < numOfNodes; ++i) {
		for (int j = 0; j < numOfNodes; ++j) {
			C_global[j][i] = 0;
		}
	}
}


void SOE::P_Clear(int numOfNodes) {
	for (int j = 0; j < numOfNodes; ++j) {
		P_global[j] = 0;
	}
	
}

void SOE::t_Print(int sizeOfT) {
	std::cout << std::endl << "T: " << std::endl;
	for (int j = 0; j < sizeOfT; ++j) {
		std::cout << t[j] << "		";
	}
}

void SOE::P_Print(int sizeOfP) {
	std::cout << std::endl << "P global: " << std::endl;
	for (int j = 0; j < sizeOfP; ++j) {
		std::cout << P_global[j] << "		";
	}
}


void SOE::H_Print(int sizeOfH) {
	std::cout << std::endl << "H global:" << std::endl;
	for (int i = 0; i < sizeOfH; ++i) {
		for (int j = 0; j < sizeOfH; ++j) {
			std::cout << H_global[j][i] << "	";
		}
		std::cout << std::endl;
	}
}

void SOE::H_Clear(int numOfNodes) {
	for (int i = 0; i < numOfNodes; ++i) {
		for (int j = 0; j < numOfNodes; ++j) {
			H_global[j][i] = 0;
		}
	}
}

void SOE::modify(Element *E,double deltaTau, int numOfElem, int numOfNodes) {
	//divide C by delta tau
	for (int wiersz = 0; wiersz < numOfNodes; ++wiersz) {
		for (int kolumna = 0; kolumna < numOfNodes; ++kolumna) {
			C_global[wiersz][kolumna] /= deltaTau;
		}
	}
	//H global + new C global
	for (int wiersz = 0; wiersz < numOfNodes; ++wiersz) {
		for (int kolumna = 0; kolumna < numOfNodes; ++kolumna) {
			H_global[wiersz][kolumna] += C_global[wiersz][kolumna];
		}
	}

	double* tmp = new double[numOfNodes];
	for (int i = 0; i < numOfNodes; ++i) {
		tmp[i] = 0;
	}

	//new C global times t0
	std::cout << std::endl << std::endl;
	for (int i = 0; i < numOfNodes; ++i) {
		for (int j = 0; j < numOfElem; ++j) {
			tmp[i] -= C_global[i][j] * t[j];		
		}
	}
	//add tmp vector to P
	for (int i = 0; i < numOfNodes; ++i) {
		P_global[i] += tmp[i];
	}
	for (int i = 0; i < numOfNodes; ++i) {
		P_global[i] *= -1;
	}
}

void SOE::genHCP_global(Element* e, int numOfElem,int numOfNodes, int numOfIntegralPoints,double alfa, double t_inf, double deltaTau) {
	H_Clear(numOfNodes);
	C_Clear(numOfNodes);
	P_Clear(numOfNodes);
	
	for (int element = 0; element < numOfElem; ++element) {
		//std::cout << "ELEMENT " << element << ": " << std::endl;
		if (numOfIntegralPoints == 2) {
			ElemUniwersal4_2point* elemUniwersal4 = new ElemUniwersal4_2point();
			elemUniwersal4->setAlfa(alfa);
			elemUniwersal4->setT_inf(t_inf);
			elemUniwersal4->genJacobians(
				e[element].getId(0)->getX(), e[element].getId(0)->getY(),
				e[element].getId(1)->getX(), e[element].getId(1)->getY(),
				e[element].getId(2)->getX(), e[element].getId(2)->getY(),
				e[element].getId(3)->getX(), e[element].getId(3)->getY()
			);
			elemUniwersal4->genH_BC(
				e[element].getId(0)->getX(), e[element].getId(0)->getY(), e[element].getId(0)->getBC(),
				e[element].getId(1)->getX(), e[element].getId(1)->getY(), e[element].getId(1)->getBC(),
				e[element].getId(2)->getX(), e[element].getId(2)->getY(), e[element].getId(2)->getBC(),
				e[element].getId(3)->getX(), e[element].getId(3)->getY(), e[element].getId(3)->getBC()
			);
			elemUniwersal4->genP(
				e[element].getId(0)->getX(), e[element].getId(0)->getY(), e[element].getId(0)->getBC(),
				e[element].getId(1)->getX(), e[element].getId(1)->getY(), e[element].getId(1)->getBC(),
				e[element].getId(2)->getX(), e[element].getId(2)->getY(), e[element].getId(2)->getBC(),
				e[element].getId(3)->getX(), e[element].getId(3)->getY(), e[element].getId(3)->getBC()
			);
			
			//elemUniwersal4->H_BCPrint();
			//elemUniwersal4->hPrint();
			//elemUniwersal4->PPrint();
			elemUniwersal4->genC(700, 7800);
			//elemUniwersal4->CPrint();
			elemUniwersal4->addH_BC2H();
			for (int row = 0; row < 4; ++row) {
				for (int column = 0; column < 4; ++column) {
					H_global[e[element].getId(row)->getGlobalId() - 1][e[element].getId(column)->getGlobalId() - 1] += elemUniwersal4->getH()[row][column];
					C_global[e[element].getId(row)->getGlobalId() - 1][e[element].getId(column)->getGlobalId() - 1] += elemUniwersal4->getC()[row][column];
				}
				P_global[e[element].getId(row)->getGlobalId() - 1] += elemUniwersal4->getP()[row]; 
			}
			//P_Print(numOfNodes);
			e[element].setH_lok(elemUniwersal4->getH(), 4, 4);
			e[element].setC_lok(elemUniwersal4->getC(), 4, 4);
			e[element].setP_lok(elemUniwersal4->getP(), 4);
			delete elemUniwersal4;
		}
		else if (numOfIntegralPoints == 3) {
			ElemUniwersal4_3point* elemUniwersal4 = new ElemUniwersal4_3point();
			elemUniwersal4->setAlfa(alfa);
			elemUniwersal4->setT_inf(t_inf);
			elemUniwersal4->genJacobians(
				e[element].getId(0)->getX(), e[element].getId(0)->getY(),
				e[element].getId(1)->getX(), e[element].getId(1)->getY(),
				e[element].getId(2)->getX(), e[element].getId(2)->getY(),
				e[element].getId(3)->getX(), e[element].getId(3)->getY()
			);
			elemUniwersal4->genH_BC(
				e[element].getId(0)->getX(), e[element].getId(0)->getY(), e[element].getId(0)->getBC(),
				e[element].getId(1)->getX(), e[element].getId(1)->getY(), e[element].getId(1)->getBC(),
				e[element].getId(2)->getX(), e[element].getId(2)->getY(), e[element].getId(2)->getBC(),
				e[element].getId(3)->getX(), e[element].getId(3)->getY(), e[element].getId(3)->getBC()
			);
			elemUniwersal4->genP(
				e[element].getId(0)->getX(), e[element].getId(0)->getY(), e[element].getId(0)->getBC(),
				e[element].getId(1)->getX(), e[element].getId(1)->getY(), e[element].getId(1)->getBC(),
				e[element].getId(2)->getX(), e[element].getId(2)->getY(), e[element].getId(2)->getBC(),
				e[element].getId(3)->getX(), e[element].getId(3)->getY(), e[element].getId(3)->getBC()
			);
			//elemUniwersal4->H_BCPrint();
			//elemUniwersal4->hPrint();
			//elemUniwersal4->PPrint();
			elemUniwersal4->genC(700, 7800);
			//elemUniwersal4->CPrint();
			elemUniwersal4->addH_BC2H();
			for (int row = 0; row < 4; ++row) {
				for (int column = 0; column < 4; ++column) {
					H_global[e[element].getId(row)->getGlobalId() - 1][e[element].getId(column)->getGlobalId() - 1] += elemUniwersal4->getH()[row][column];
					C_global[e[element].getId(row)->getGlobalId() - 1][e[element].getId(column)->getGlobalId() - 1] += elemUniwersal4->getC()[row][column];
				}
				P_global[e[element].getId(row)->getGlobalId() - 1] += elemUniwersal4->getP()[row];
			}

			e[element].setH_lok(elemUniwersal4->getH(), 4, 4);
			e[element].setC_lok(elemUniwersal4->getC(), 4, 4);
			e[element].setP_lok(elemUniwersal4->getP(), 4);
			delete elemUniwersal4;
		}
		else {
			std::cout << "Wrong integral number" << std::endl;
			return;
		}

	}
	//std::cout << "Przed modify" << std::endl;
	//H_Print(numOfNodes);
	//C_Print(numOfNodes);
	//P_Print(numOfNodes);

	modify(e, deltaTau, numOfElem, numOfNodes);
	
	//std::cout << "Po modify" << std::endl;
	//H_Print(numOfNodes);
	//C_Print(numOfNodes);
	//P_Print(numOfNodes);

	
	
}


SOE::SOE() {
	H_global = new double*[1];
	H_global[0] = new double[1];
	C_global = new double*[1];
	C_global[0] = new double[1];
	P_global = new double[1];
	t = new double[1];
}

SOE::SOE(int sizeOfH) {
	H_global = new double *[sizeOfH];
	for (int i = 0; i < sizeOfH; ++i) {
		H_global[i] = new double[sizeOfH];
	}

	C_global = new double *[sizeOfH];
	for (int i = 0; i < sizeOfH; ++i) {
		C_global[i] = new double[sizeOfH];
	}
	P_global = new double[sizeOfH];
	t = new double[sizeOfH];
}

double SOE::getElemFromH(int row, int column) {
	return H_global[row][column];
}

void SOE::setElemFromH(int row, int column, double value) {
	H_global[row][column] = value;
}

void SOE::addElemToH(int row, int column, double value) {
	H_global[row][column] += value;
}

GlobalData::GlobalData() {
	std::fstream file;
	std::string text;


	file.open("INPUT.txt", std::ios::in);
	//setting GlobalData params
	if (file.good()) {
		//H
		getline(file, text);
		H = atof(text.c_str());
		//W
		getline(file, text);
		W = atof(text.c_str());
		//nH
		getline(file, text);
		nH = atoi(text.c_str());
		//nW
		getline(file, text);
		nW = atoi(text.c_str());
		//nN
		nN = nW * nH;
		//nE
		nE = (nH - 1) * (nW - 1);
	}
	else std::cout << "file_read ERROR" << std::endl;

	file.close();
}

void GlobalData::setBC(Element* E) {
	//left side
	for (int i = 0; i <= nH - 2; ++i) {
		E[i].getId(0)->setBC(1);
		E[i].getId(3)->setBC(1);
	}
	
	//E[nH - 2].getId(3)->setBC(1);		//top right corner
	//bottom side
	for (int i = 0; i < nE; i += nH - 1) {
		
		E[i].getId(0)->setBC(1);
		E[i].getId(1)->setBC(1);
		
		
		/*//bottom
		if ((i + 1) % nH == 0) {
			E[i].getId(0)->setBC(1);
		}
		//top
		if (i % nH == 1) {
			E[i].getId(4)->setBC(1);
		}*/
	}
	//top side
	for (int i = nH - 1; i <= nE; i += nH - 1) {
		E[i - 1].getId(3)->setBC(1);
		E[i - 1].getId(2)->setBC(1);
	}

	//right side
	for (int i = nE - nH + 1; i <= nE - 1; ++i) {
		E[i].getId(1)->setBC(1);
		E[i].getId(2)->setBC(1);
	}
	
}


void GlobalData::genArrays(Element* E, Node* N,double t0) {
	// deltaH & deltaW distance between two nodes
	double deltaH = this->getH() / (this->getNH() - 1);
	double deltaW = this->getW() / (this->getNW() - 1);

	//gen nodes
	double x = 0, y = 0;
	int nodeId = 1;
	for (int i = 0; i < this->getNN(); ++i) {
		if (y <= this->getH()) {
			N[i].setX(x);
			N[i].setY(y);
			N[i].setGlobalId(nodeId);
			N[i].setT(t0);
			nodeId++;
			y += deltaH;
		}
		else {
			x += deltaW;
			y = 0;
			N[i].setX(x);
			N[i].setY(y);
			N[i].setGlobalId(nodeId);
			N[i].setT(t0);
			nodeId++;
			y += deltaH;
		}
	}
	//gen elements
	int tmp = 1;
	int adder = 0;
	int counter = 0;

	for (int i = 0; i < this->getNE() - counter; ++i) {
		
		if (tmp > this->getNH() - 1) {
			tmp = 1;
			counter--;
			continue;
		}

		E[i + counter].setId(N[i + adder], 0);
		E[i + counter].setId(N[i + 1 + adder], 3);
		E[i + counter].setId(N[i + this->getNH() + adder], 1);
		E[i + counter].setId(N[i + 1 + this->getNH() + adder], 2);

		/*
		std::cout << "Element " << i << std::endl;
		std::cout << "1. id = " << N[i].getGlobalId() << "	x = " << N[i].getX() << "		y = " << N[i].getY() << std::endl;
		std::cout << "2. id = " << N[i + 1].getGlobalId() << "	x = " << N[i + 1].getX() << "		y = " << N[i + 1].getY() << std::endl;
		std::cout << "3. id = " << N[i + this->getNH()].getGlobalId() << "	x = " << N[i + this->getNH()].getX() << "		y = " << N[i + this->getNH()].getY() << std::endl;
		std::cout << "4. id = " << N[i + this->getNH() + 1].getGlobalId() << "	x = " << N[i + this->getNH() + 1].getX() << "		y = " << N[i + this->getNH() + 1].getY() << std::endl;
		*/
		
		tmp++;
	}

	setBC(E);
}

#endif // !NET_GEN
