#include "pch.h"
#include "Koef.h"
//#include "Koef.h"

//#include <math.h>
//#include "Constants.h"

double DampingKoef::linInterp(double x1, double x2, double y1, double y2, double x) {
	double y = (y2 * (x - x1) + y1 * (x2 - x)) / (x2 - x1);
	return(y);

}

DampingKoef::DampingKoef()
{
}

double DampingKoef::MrFr(double mach, double alpha) {
	double alphaToDeg = abs(alpha) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMrFr[0][j + 1]; j++);
	double MrFr = 0;
	if (mach < 0.6) {
		MrFr = linInterp(ArMrFr[0][j], ArMrFr[0][j+1], ArMrFr[1][j], ArMrFr[1][j+1], alphaToDeg);
	} else {
		if (mach > 2) {
			MrFr = linInterp(ArMrFr[0][j], ArMrFr[0][j+1], ArMrFr[4][j], ArMrFr[4][j+1], alphaToDeg);
		} else {
			int i;
			for (i = 2; (mach > ArMrFr[i][0]); i++);
			double Mr1 = linInterp(ArMrFr[i - 1][0], ArMrFr[i][0], ArMrFr[i - 1][j], ArMrFr[i][j], mach);
			double Mr2 = linInterp(ArMrFr[i - 1][0], ArMrFr[i][0], ArMrFr[i - 1][j+1], ArMrFr[i][j+1], mach);
			MrFr = linInterp(ArMrFr[0][j - 1], ArMrFr[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MrFr;
}
double DampingKoef::MpFr(double mach, double alpha) {
	double alphaToDeg = abs(alpha) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMpFr[0][j + 1]; j++);
	double MpFr = 0;
	if (mach < 0.6) {
		MpFr = linInterp(ArMpFr[0][j], ArMpFr[0][j + 1], ArMpFr[1][j], ArMpFr[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MpFr = linInterp(ArMpFr[0][j], ArMpFr[0][j + 1], ArMpFr[4][j], ArMpFr[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMpFr[i][0]); i++);
			double Mr1 = linInterp(ArMpFr[i - 1][0], ArMpFr[i][0], ArMpFr[i - 1][j], ArMpFr[i][j], mach);
			double Mr2 = linInterp(ArMpFr[i - 1][0], ArMpFr[i][0], ArMpFr[i - 1][j + 1], ArMpFr[i][j + 1], mach);
			MpFr = linInterp(ArMpFr[0][j - 1], ArMpFr[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MpFr;
}
double DampingKoef::MyFr(double mach, double alpha) {
	double alphaToDeg = abs(alpha) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMyFr[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(ArMyFr[0][j], ArMyFr[0][j + 1], ArMyFr[1][j], ArMyFr[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MyFr = linInterp(ArMyFr[0][j], ArMyFr[0][j + 1], ArMyFr[4][j], ArMyFr[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMyFr[i][0]); i++);
			double Mr1 = linInterp(ArMyFr[i - 1][0], ArMyFr[i][0], ArMyFr[i - 1][j], ArMyFr[i][j], mach);
			double Mr2 = linInterp(ArMyFr[i - 1][0], ArMyFr[i][0], ArMyFr[i - 1][j + 1], ArMyFr[i][j + 1], mach);
			MyFr = linInterp(ArMyFr[0][j - 1], ArMyFr[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}


double DampingKoef::MrFp(double mach, double alpha) {
	double alphaToDeg = abs(alpha) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMrFp[0][j + 1]; j++);
	double MrFr = 0;
	if (mach < 0.6) {
		MrFr = linInterp(ArMrFp[0][j], ArMrFp[0][j + 1], ArMrFp[1][j], ArMrFp[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MrFr = linInterp(ArMrFp[0][j], ArMrFp[0][j + 1], ArMrFp[4][j], ArMrFp[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMrFp[i][0]); i++);
			double Mr1 = linInterp(ArMrFp[i - 1][0], ArMrFp[i][0], ArMrFp[i - 1][j], ArMrFp[i][j], mach);
			double Mr2 = linInterp(ArMrFp[i - 1][0], ArMrFp[i][0], ArMrFp[i - 1][j + 1], ArMrFp[i][j + 1], mach);
			MrFr = linInterp(ArMrFp[0][j - 1], ArMrFp[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MrFr;
}
double DampingKoef::MpFp(double mach, double alpha) {
	double alphaToDeg = abs(alpha) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMpFp[0][j + 1]; j++);
	double MpFr = 0;
	if (mach < 0.6) {
		MpFr = linInterp(ArMpFp[0][j], ArMpFp[0][j + 1], ArMpFp[1][j], ArMpFp[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MpFr = linInterp(ArMpFp[0][j], ArMpFp[0][j + 1], ArMpFp[4][j], ArMpFp[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMpFp[i][0]); i++);
			double Mr1 = linInterp(ArMpFp[i - 1][0], ArMpFp[i][0], ArMpFp[i - 1][j], ArMpFp[i][j], mach);
			double Mr2 = linInterp(ArMpFp[i - 1][0], ArMpFp[i][0], ArMpFp[i - 1][j + 1], ArMpFp[i][j + 1], mach);
			MpFr = linInterp(ArMpFp[0][j - 1], ArMpFp[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MpFr;
}
double DampingKoef::MyFp(double mach, double alpha) {
	double alphaToDeg = abs(alpha) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMyFp[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(ArMyFp[0][j], ArMyFp[0][j + 1], ArMyFp[1][j], ArMyFp[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MyFr = linInterp(ArMyFp[0][j], ArMyFp[0][j + 1], ArMyFp[4][j], ArMyFp[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMyFp[i][0]); i++);
			double Mr1 = linInterp(ArMyFp[i - 1][0], ArMyFp[i][0], ArMyFp[i - 1][j], ArMyFp[i][j], mach);
			double Mr2 = linInterp(ArMyFp[i - 1][0], ArMyFp[i][0], ArMyFp[i - 1][j + 1], ArMyFp[i][j + 1], mach);
			MyFr = linInterp(ArMyFp[0][j - 1], ArMyFp[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}

double DampingKoef::MrFy(double mach, double alpha) {
	double alphaToDeg = abs(alpha) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMrFy[0][j + 1]; j++);
	double MrFr = 0;
	if (mach < 0.6) {
		MrFr = linInterp(ArMrFy[0][j], ArMrFy[0][j + 1], ArMrFy[1][j], ArMrFy[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MrFr = linInterp(ArMrFy[0][j], ArMrFy[0][j + 1], ArMrFy[4][j], ArMrFy[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMrFy[i][0]); i++);
			double Mr1 = linInterp(ArMrFy[i - 1][0], ArMrFy[i][0], ArMrFy[i - 1][j], ArMrFy[i][j], mach);
			double Mr2 = linInterp(ArMrFy[i - 1][0], ArMrFy[i][0], ArMrFy[i - 1][j + 1], ArMrFy[i][j + 1], mach);
			MrFr = linInterp(ArMrFy[0][j - 1], ArMrFy[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MrFr;
}
double DampingKoef::MpFy(double mach, double alpha) {
	double alphaToDeg = abs(alpha) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMpFy[0][j + 1]; j++);
	double MpFr = 0;
	if (mach < 0.6) {
		MpFr = linInterp(ArMpFy[0][j], ArMpFy[0][j + 1], ArMpFy[1][j], ArMpFy[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MpFr = linInterp(ArMpFy[0][j], ArMpFy[0][j + 1], ArMpFy[4][j], ArMpFy[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMpFy[i][0]); i++);
			double Mr1 = linInterp(ArMpFy[i - 1][0], ArMpFy[i][0], ArMpFy[i - 1][j], ArMpFy[i][j], mach);
			double Mr2 = linInterp(ArMpFy[i - 1][0], ArMpFy[i][0], ArMpFy[i - 1][j + 1], ArMpFy[i][j + 1], mach);
			MpFr = linInterp(ArMpFy[0][j - 1], ArMpFy[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MpFr;
}
double DampingKoef::MyFy(double mach, double alpha) {
	double alphaToDeg = abs(alpha) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMyFy[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(ArMyFy[0][j], ArMyFy[0][j + 1], ArMyFy[1][j], ArMyFy[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MyFr = linInterp(ArMyFy[0][j], ArMyFy[0][j + 1], ArMyFy[4][j], ArMyFy[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMyFy[i][0]); i++);
			double Mr1 = linInterp(ArMyFy[i - 1][0], ArMyFy[i][0], ArMyFy[i - 1][j], ArMyFy[i][j], mach);
			double Mr2 = linInterp(ArMyFy[i - 1][0], ArMyFy[i][0], ArMyFy[i - 1][j + 1], ArMyFy[i][j + 1], mach);
			MyFr = linInterp(ArMyFy[0][j - 1], ArMyFy[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}

double DampingKoef::MrB(double mach, double beta) {
	double alphaToDeg = abs(beta) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMrFb[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(ArMrFb[0][j], ArMrFb[0][j + 1], ArMrFb[1][j], ArMrFb[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 1.2) {
			MyFr = linInterp(ArMrFb[0][j], ArMrFb[0][j + 1], ArMrFb[3][j], ArMrFb[3][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMrFb[i][0]); i++);
			double Mr1 = linInterp(ArMrFb[i - 1][0], ArMrFb[i][0], ArMrFb[i - 1][j], ArMrFb[i][j], mach);
			double Mr2 = linInterp(ArMrFb[i - 1][0], ArMrFb[i][0], ArMrFb[i - 1][j + 1], ArMrFb[i][j + 1], mach);
			MyFr = linInterp(ArMrFb[0][j - 1], ArMrFb[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}
double DampingKoef::MpB(double mach, double beta) {
	double alphaToDeg = abs(beta) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMpFb[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(ArMpFb[0][j], ArMpFb[0][j + 1], ArMpFb[1][j], ArMpFb[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 1.2) {
			MyFr = linInterp(ArMpFb[0][j], ArMpFb[0][j + 1], ArMpFb[3][j], ArMpFb[3][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMpFb[i][0]); i++);
			double Mr1 = linInterp(ArMpFb[i - 1][0], ArMpFb[i][0], ArMpFb[i - 1][j], ArMpFb[i][j], mach);
			double Mr2 = linInterp(ArMpFb[i - 1][0], ArMpFb[i][0], ArMpFb[i - 1][j + 1], ArMpFb[i][j + 1], mach);
			MyFr = linInterp(ArMpFb[0][j - 1], ArMpFb[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}
double DampingKoef::MyB(double mach, double beta) {
	double alphaToDeg = abs(beta) * toDeg;
	int j;
	for (j = 1; alphaToDeg > ArMyFb[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(ArMyFb[0][j], ArMyFb[0][j + 1], ArMyFb[1][j], ArMyFb[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 1.2) {
			MyFr = linInterp(ArMyFb[0][j], ArMyFb[0][j + 1], ArMyFb[3][j], ArMyFb[3][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > ArMyFb[i][0]); i++);
			double Mr1 = linInterp(ArMyFb[i - 1][0], ArMyFb[i][0], ArMyFb[i - 1][j], ArMyFb[i][j], mach);
			double Mr2 = linInterp(ArMyFb[i - 1][0], ArMyFb[i][0], ArMyFb[i - 1][j + 1], ArMyFb[i][j + 1], mach);
			MyFr = linInterp(ArMyFb[0][j - 1], ArMyFb[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}

DampingKoef::~DampingKoef()
{
}
