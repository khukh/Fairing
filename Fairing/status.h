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

	void addErotation(double lon, double azim);

	double pressure;

	const double R_EARTH_G = (6371E3);
	const double PI0 = (398600.4418E9);


	Atmosphere GOST4401;
	Vect<3> Fg; //проекция в-ра сил на старт СК
	Vect<3> ForcePr;
	//std::vector <double> ForcePrG;  //проекция 
	//std::vector <double> ForceG;  //вектор силы притяжения
	Vect<3> Torque, v;  //момент силы

	rotation Rot;


	virtual void nonIntegr();  //пересчет неинтегрируемых параметров
	virtual Vect<15> rightPart();	//значения производных	

	double getH();
	double vv;
	double getV();
	void setStageParam(double mDry, double mFuel, double mp, double sm, double sa, double l, double p, double ix, double iy, double iz);


	void addV(double Vx, double Vy, double vz);


	void printParam( std::ofstream &fout);	//вывод параметров
	void setParam(Vect<15> b);
	Vect<15> getParam();



	status& operator=(const status& right);
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

	double alpha;
	double betta;
	//double S_M = PI *D_M*D_M/4;
	double m;//масса


	double g;

	double density, q, mach, vFullsq;
	double P;
	double M_P;

	//параметры ракеты
	double M_DRY, M_FUEL;
	double S_M, S_A;
	double L;
	double I_X;
	double I_Y;
	double I_Z;

	double azimut, lon, lat;
	double Ve;

	double mzAl, h, cx, cy, cz, mzBet;




};

