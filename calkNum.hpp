#ifndef CALK_NUM
#define CALK_NUM
#include<math.h>

typedef double(*function)(double, double);
/*
double function(double x, double y) {
	return x + y;
}
*/
double calka2DStandarized2Points(function function) {
	double pierw3 = sqrt(3);

	double pcIx[4] = { (-1 / pierw3), (1 / pierw3), (1 / pierw3), (-1 / pierw3) };
	double pcIy[4] = { (-1 / pierw3), (-1 / pierw3), (1 / pierw3),(1 / pierw3) };
	double retValue = 0;

	for (int i = 0; i < 4; ++i) {
		retValue += function(pcIx[i], pcIy[i]);
	}
	return retValue;
}

double calka2DStandarized3Points(function function) {
	double pierw3_5 = 0.7745966692414834;
	double pc[3] = { -pierw3_5, 0, pierw3_5 };
	double w[3] = { 0.5555555555555555 , 0.8888888888888888, 0.5555555555555555 };

	double retValue = 0;

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			retValue += (function(pc[j], pc[i]) * w[i] * w[j]);
		}
	}
	return retValue;
}

#endif // !CALK_NUM
