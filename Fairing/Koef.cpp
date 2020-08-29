#include "pch.h"
#include "Koef.h"
//#include "Koef.h"

//#include <math.h>
//#include "Constants.h"

double DampingKoef::linInterp(double x1, double x2, double y1, double y2, double x) {
	double y = (y2 * (x - x1) + y1 * (x2 - x)) / (x2 - x1);
	return(y);
}


double DampingKoef::MrFr(double mach, double alpha) {
	double Mr[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, -0.425, -0.45, 0, -0.75, -0.1866, 3.8134, -0.2134, 1.92, -0.425},
		{0.98, -0.4, -0.7, -0.775, -0.5, -0.2666, 0.75, 0.025, -4.75, -0.4},
		{1.2, -0.4, -0.8, 0.6, -0.425, -0.2666, -0.1066, 0, -0.65, -0.4},
		{2, -0.4, 0.025, -0.5, -0.6, -0.2934, -1.0934, 0.8266, -0.8, -0.4}
	};
	double alphaToDeg = alpha * toDeg;
	int j;
	for (j = 1; alphaToDeg > Mr[0][j + 1]; j++);
	double MrFr = 0;
	if (mach < 0.6) {
		MrFr = linInterp(Mr[0][j], Mr[0][j+1], Mr[1][j], Mr[1][j+1], alphaToDeg);
	} else {
		if (mach > 2) {
			MrFr = linInterp(Mr[0][j], Mr[0][j+1], Mr[4][j], Mr[4][j+1], alphaToDeg);
		} else {
			int i;
			for (i = 2; (mach > Mr[i][0]); i++);
			double Mr1 = linInterp(Mr[i - 1][0], Mr[i][0], Mr[i - 1][j], Mr[i][j], mach);
			double Mr2 = linInterp(Mr[i - 1][0], Mr[i][0], Mr[i - 1][j+1], Mr[i][j+1], mach);
			MrFr = linInterp(Mr[0][j - 1], Mr[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MrFr;
}
double DampingKoef::MpFr(double mach, double alpha) {
	double Mp[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, 0, 0, -0.3692, -0.1846, -0.0308, -0.0308, 1.5692, 0.3076, 0},
		{0.98, -0.0308, 0.0616, 1.2308, 0.3384, -0.0308, -0.0308, -0.123, 0.3692, -0.0308},
		{1.2, 0, 0.2154, -1.0154, -1.877, 0, -0.0616, -4.9538, -0.277, 0},
		{2, 0, -0.3384, -0.8, -2.4616, -0.0616, -0.0308, 4.2154, 0.523, 0}
	};
	double alphaToDeg = alpha * toDeg;
	int j;
	for (j = 1; alphaToDeg > Mp[0][j + 1]; j++);
	double MpFr = 0;
	if (mach < 0.6) {
		MpFr = linInterp(Mp[0][j], Mp[0][j + 1], Mp[1][j], Mp[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MpFr = linInterp(Mp[0][j], Mp[0][j + 1], Mp[4][j], Mp[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > Mp[i][0]); i++);
			double Mr1 = linInterp(Mp[i - 1][0], Mp[i][0], Mp[i - 1][j], Mp[i][j], mach);
			double Mr2 = linInterp(Mp[i - 1][0], Mp[i][0], Mp[i - 1][j + 1], Mp[i][j + 1], mach);
			MpFr = linInterp(Mp[0][j - 1], Mp[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MpFr;
}
double DampingKoef::MyFr(double mach, double alpha) {
	double My[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, -0.338, -1.446, -0.892, 0.431, 0.308, 3.2, 1.015, -4.954, -0.246},
		{0.98, -0.154, -0.892, -0.031, 1.292, 0.462, 2.615, 2.215, -1.446, -0.123},
		{1.2, -0.277, -0.646, 0.431, -0.092, 0.646, 1.323, 0.154, -2.615, -0.215},
		{2, -0.277, 0.031, 0.4, 0.985, 0.277, 0.154, 1.015, -1.477, -0.369}
	};
	double alphaToDeg = alpha * toDeg;
	int j;
	for (j = 1; alphaToDeg > My[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(My[0][j], My[0][j + 1], My[1][j], My[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MyFr = linInterp(My[0][j], My[0][j + 1], My[4][j], My[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > My[i][0]); i++);
			double Mr1 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j], My[i][j], mach);
			double Mr2 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j + 1], My[i][j + 1], mach);
			MyFr = linInterp(My[0][j - 1], My[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}


double DampingKoef::MrFp(double mach, double alpha) {
	double Mr[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, 0, 0, 0.1412, -0.2625, -0.0187, -0.075, -0.6, 0.2471, 0.0},
		{0.98, 0.0, -0.0187, 0.4059, -0.4235, 0, -0.3, -2.8125, 0.0529, 0},
		{1.2, 0.0, -0.1687, -1.0412, -0.3, 0, 0.0529, -0.9882, 0.7412, 0},
		{2, 0.0, 0.5471, -0.075, -0.1875, 0, 0.1412, -3.0176, 0.4588, 0}
	};
	double alphaToDeg = alpha * toDeg;
	int j;
	for (j = 1; alphaToDeg > Mr[0][j + 1]; j++);
	double MrFr = 0;
	if (mach < 0.6) {
		MrFr = linInterp(Mr[0][j], Mr[0][j + 1], Mr[1][j], Mr[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MrFr = linInterp(Mr[0][j], Mr[0][j + 1], Mr[4][j], Mr[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > Mr[i][0]); i++);
			double Mr1 = linInterp(Mr[i - 1][0], Mr[i][0], Mr[i - 1][j], Mr[i][j], mach);
			double Mr2 = linInterp(Mr[i - 1][0], Mr[i][0], Mr[i - 1][j + 1], Mr[i][j + 1], mach);
			MrFr = linInterp(Mr[0][j - 1], Mr[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MrFr;
}
double DampingKoef::MpFp(double mach, double alpha) {
	double Mp[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, -5.187, -3.882, -0.562, -2.625, -6.294, -8.75, -1.625, -9.812, -5.187},
		{0.98, -9, -5.25, 1.438, -1.125, -10, -12.375, 1.938, -10.706, -9},
		{1.2, -10.706, -4.937, -6.176, -2.625, -10.471, -9.875, -4.187, -10.824, -10.706},
		{2, -8.25, -5.812, -8.5, -8.937, -5.812, -14, -13.353, -11.125, -8.25}
	};
	double alphaToDeg = alpha * toDeg;
	int j;
	for (j = 1; alphaToDeg > Mp[0][j + 1]; j++);
	double MpFr = 0;
	if (mach < 0.6) {
		MpFr = linInterp(Mp[0][j], Mp[0][j + 1], Mp[1][j], Mp[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MpFr = linInterp(Mp[0][j], Mp[0][j + 1], Mp[4][j], Mp[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > Mp[i][0]); i++);
			double Mr1 = linInterp(Mp[i - 1][0], Mp[i][0], Mp[i - 1][j], Mp[i][j], mach);
			double Mr2 = linInterp(Mp[i - 1][0], Mp[i][0], Mp[i - 1][j + 1], Mp[i][j + 1], mach);
			MpFr = linInterp(Mp[0][j - 1], Mp[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MpFr;
}
double DampingKoef::MyFp(double mach, double alpha) {
	double My[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, -0.0, -0.0159, -0.1263, -0.4737, -0.0474, -0.0474, -0.4896, 0.3474, 0},
		{0.98, -0.0, -0.0789, 0.2211, -1.1367, -0.0159, -0.6159, -3.6159, -0.0474, 0.0},
		{1.2, -0.0, -0.2211, 0.1578, 1.3896, 0.0474, 0, -0.1737, -0.2211, -0.0},
		{2, -0.0, 0.3315, -0.0315, 0.8685, 0.0315, -0.0474, -1.1526, 1.5, -0.0}

	};
	double alphaToDeg = alpha * toDeg;
	int j;
	for (j = 1; alphaToDeg > My[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(My[0][j], My[0][j + 1], My[1][j], My[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MyFr = linInterp(My[0][j], My[0][j + 1], My[4][j], My[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > My[i][0]); i++);
			double Mr1 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j], My[i][j], mach);
			double Mr2 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j + 1], My[i][j + 1], mach);
			MyFr = linInterp(My[0][j - 1], My[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}

double DampingKoef::MrFy(double mach, double alpha) {
	double Mr[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, 1.686, -0.232, 0.371, -0.087, -1.536, -6.265, -0.319, 4.629, 1.686},
		{0.98, 1.257, -1.217, -0.261, 1.171, -0.783, -2.171, -0.29, -3.6, 1.257},
		{1.2, 1.229, -0.377, -0.029, 0.571, -0.464, -1.275, -0.435, 0.6, 1.229},
		{2, 1.486, 0.114, -0.29, 0.143, -0.58, -0.638, 0.457, 0.829, 1.486}
	};
	double alphaToDeg = alpha * toDeg;
	int j;
	for (j = 1; alphaToDeg > Mr[0][j + 1]; j++);
	double MrFr = 0;
	if (mach < 0.6) {
		MrFr = linInterp(Mr[0][j], Mr[0][j + 1], Mr[1][j], Mr[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MrFr = linInterp(Mr[0][j], Mr[0][j + 1], Mr[4][j], Mr[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > Mr[i][0]); i++);
			double Mr1 = linInterp(Mr[i - 1][0], Mr[i][0], Mr[i - 1][j], Mr[i][j], mach);
			double Mr2 = linInterp(Mr[i - 1][0], Mr[i][0], Mr[i - 1][j + 1], Mr[i][j + 1], mach);
			MrFr = linInterp(Mr[0][j - 1], Mr[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MrFr;
}
double DampingKoef::MpFy(double mach, double alpha) {
	double Mp[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, -0.516, -0.548, -0.548, -0.6915, -0.468, -0.484, -0.452, -0.3405, -0.516},
		{0.98, -0.5, -0.532, -0.564, 0.8885, -0.516, -0.165, -0.532, -0.2285, -0.5},
		{1.2, -0.5, -0.819, -0.149, -2.516, -0.532, -0.436, -0.4045, -0.532, -0.5},
		{2, -0.516, -0.947, -0.787, -2.8225, -0.516, -0.452, 2.3245, 0.1385, -0.516}
	};
	double alphaToDeg = alpha * toDeg;
	int j;
	for (j = 1; alphaToDeg > Mp[0][j + 1]; j++);
	double MpFr = 0;
	if (mach < 0.6) {
		MpFr = linInterp(Mp[0][j], Mp[0][j + 1], Mp[1][j], Mp[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MpFr = linInterp(Mp[0][j], Mp[0][j + 1], Mp[4][j], Mp[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > Mp[i][0]); i++);
			double Mr1 = linInterp(Mp[i - 1][0], Mp[i][0], Mp[i - 1][j], Mp[i][j], mach);
			double Mr2 = linInterp(Mp[i - 1][0], Mp[i][0], Mp[i - 1][j + 1], Mp[i][j + 1], mach);
			MpFr = linInterp(Mp[0][j - 1], Mp[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MpFr;
}
double DampingKoef::MyFy(double mach, double alpha) {
	double My[5][10] = {
		{0, -180, -135, -90, -45, 0, 45, 90, 135, 180},
		{0.6, -3.0937, -2.6512, -6.7442, -5.8594, -3.2578, -2.186, -0.1395, -3.2344, -3.0937},
		{0.98, -3.5859, -3.4922, -4.5234, -4.5234, -4.5234, -3.4687, 0.4385, -1.3256, -3.5859},
		{1.2, -3.9375, -3.5859, -4.125, -5.1094, -4.5937, -3.0703, 0.2538, -2.7442, -3.9375},
		{2, -4.2891, -2.3488, -3.6562, -3.1172, -3.0937, -2.7442, -1.4186, -2.093, -4.2891}


	};
	double alphaToDeg = alpha * toDeg;
	int j;
	for (j = 1; alphaToDeg > My[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(My[0][j], My[0][j + 1], My[1][j], My[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 2) {
			MyFr = linInterp(My[0][j], My[0][j + 1], My[4][j], My[4][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > My[i][0]); i++);
			double Mr1 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j], My[i][j], mach);
			double Mr2 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j + 1], My[i][j + 1], mach);
			MyFr = linInterp(My[0][j - 1], My[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}

double DampingKoef::MrB(double mach, double beta) {
	double My[4][8] = {
		{0, -90, -60, -30, 0, 30, 60, 90},
		{0.6, 1.333, 2.604, 0.125, -0.208, 0.125, 2.625, 1.354},
		{0.98, -0.354, 0.542, -0.354, -0.292, -0.333, 0.542, -0.396},
		{1.2, -0.792, -0.375, -0.729, -0.375, -0.75, -0.396, -0.792}


	};
	double alphaToDeg = beta * toDeg;
	int j;
	for (j = 1; alphaToDeg > My[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(My[0][j], My[0][j + 1], My[1][j], My[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 1.2) {
			MyFr = linInterp(My[0][j], My[0][j + 1], My[3][j], My[3][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > My[i][0]); i++);
			double Mr1 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j], My[i][j], mach);
			double Mr2 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j + 1], My[i][j + 1], mach);
			MyFr = linInterp(My[0][j - 1], My[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}
double DampingKoef::MpB(double mach, double beta) {
	double My[4][8] = {
		{0, -90, -60, -30, 0, 30, 60, 90},
		{0.6, -1.84, 1.85, -2.519, -3.897, -2.519, 1.85, -1.783},
		{0.98, -2.745, -0.028, -4.935, -8.355, -4.935, 0, -2.774},
		{1.2, -3.28, -2.972, -4.991, -8.888, -4.935, -2.943, -3.308}


	};
	double alphaToDeg = beta * toDeg;
	int j;
	for (j = 1; alphaToDeg > My[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(My[0][j], My[0][j + 1], My[1][j], My[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 1.2) {
			MyFr = linInterp(My[0][j], My[0][j + 1], My[3][j], My[3][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > My[i][0]); i++);
			double Mr1 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j], My[i][j], mach);
			double Mr2 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j + 1], My[i][j + 1], mach);
			MyFr = linInterp(My[0][j - 1], My[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}
double DampingKoef::MyB(double mach, double beta) {
	double My[4][8] = {
		{0, -90, -60, -30, 0, 30, 60, 90},
		{0.6, 2.477, 0.209, -3.244, -3.244, -3.349, 0.174, 2.512},
		{0.98, 0.663, -3.593, -4.744, -4.605, -4.779, -3.628, 0.663},
		{1.2, 0.628, -2.442, -4.535, -4.5, -4.535, -2.477, 0.523},


	};
	double alphaToDeg = beta * toDeg;
	int j;
	for (j = 1; alphaToDeg > My[0][j + 1]; j++);
	double MyFr = 0;
	if (mach < 0.6) {
		MyFr = linInterp(My[0][j], My[0][j + 1], My[1][j], My[1][j + 1], alphaToDeg);
	}
	else {
		if (mach > 1.2) {
			MyFr = linInterp(My[0][j], My[0][j + 1], My[3][j], My[3][j + 1], alphaToDeg);
		}
		else {
			int i;
			for (i = 2; (mach > My[i][0]); i++);
			double Mr1 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j], My[i][j], mach);
			double Mr2 = linInterp(My[i - 1][0], My[i][0], My[i - 1][j + 1], My[i][j + 1], mach);
			MyFr = linInterp(My[0][j - 1], My[0][j], Mr1, Mr2, alphaToDeg);
		}
	}
	//std::cout << i << '\t' << j << '\t' << cx << '\n';
	return MyFr;
}

DampingKoef::~DampingKoef()
{
}

DampingKoef::DampingKoef()
{
}
