#pragma once
class Atmosphere {
public:

	const double R = 287.05287;
	const double GST = 9.80665;
	const double PST = 101325.0;
	const double TST = 288.15;
	const double BS = 1.458E-6;
	const double S = 110.4;
	


	Atmosphere();

	double pFunc(double h);
	double tFunc(double h);
	double roFunc(double h);
	double aFunc(double h);
	double muFunc(double h);
	double nuFunc(double h);
	
	~Atmosphere();
	//double tfunc(double h);
private:
	void fillArrPandT();
	double hFromgH(double gH);
	double gHFh(double h);
	double betaKoef(int i);
	int index(double geoPotH);
	int i;

	//значения геопотенциальной высоты - первая строчка массива в м
	//значения коэффициента бета=dT/dH - вторая строка массива
	double Beta[2][12] = {
		{ -2000.0, 0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 85000.0, 94.0E+3, 102.45E+3, 117.777E+3 },
		{ -6.5E-3, -6.5E-3, 0.0, 1.0E-3, 2.8E-3, 0.0, -2.8E-3, -2.0E-3, 0.0, 3E-3, 11E-3 }
	};
	double gPotHStand[12] = { -2.0E+3, 0.0, 11.0E+3, 20.0E+3, 32.0E+3, 47.0E+3, 51.0E+3, 71.0E+3, 85.0E+3, 94.0E+3, 102.45E+3, 117.777E+3 };

	double pGostStand[12];
	double temperature[12];
};

