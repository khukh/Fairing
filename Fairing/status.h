#pragma once

#include "pch.h"

#include "iostream" 
#include <fstream>
#include "rotation.h"
#include "Atmosphere.h"
#include "KoefCalc.h"



class status {
public:
	status(std::vector <double> &a);
	~status();
	//methods
	void addErotation(double lon, double azim);
	void setStageParam(double mDry, double mFuel, double mp, double sm, double sa, double l, double p, double ix, double iy, double iz);
	virtual void nonIntegr();  //пересчет неинтегрируемых параметров
	void addV(double Vx, double Vy, double vz);
	void printParam(std::ofstream &fout);	//вывод параметров
	void setParam(Vect<15> b);
	virtual Vect<15> rightPart();	//значения производных	
	Vect<15> getParam();
	status& operator=(const status& right);
	double getH();
	double getV();




	
protected:

	Vect<15> parametr;
	/*	parametr [0] = Vx
		parametr [1] = Vy
		parametr [2] = Vz
		parametr [3] = x
		parametr [4] = y
		parametr [5] = z
		parametr [6] = Wx
		parametr [7] = Wy
		parametr [8] = Wz
		parametr [9] = ro
		parametr [10] = l
		parametr [11] = mu
		parametr [12] = nu
		parametr [13] = t
		parametr [14] = m

		*/

	rotation Rot;
	double azimut, lon, lat;
	double h;
	double Ve;
	Vect<3> v;
	double vv, vFullsq;


	Atmosphere GOST4401;
	double density, pressure;
	double q, mach;

	double alpha, betta;
	double fi1, al1;

	double cx, cy, cz;
	double mzwz, mzAl, myBet, mx;

	Vect<3> ForcePr;
	Vect<3> Fg; //проекция в-ра сил на старт СК

	const double R_EARTH_G = (6371E3);
	const double PI0 = (398600.4418E9);
	double g;
	Vect<3> Torque;  //момент силы
	double S_M, S_A;
	double L;
	double I_X;
	double I_Y;
	double I_Z;



	double m;//масса


	


	double P;
	double M_P;

	//параметры ракеты
	double M_DRY, M_FUEL;


	double d1, a1, p1;

};

