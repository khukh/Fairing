
#include "pch.h"

#define _USE_MATH_DEFINES
#include "status.h"

#include "cmath"
#include <fstream>


status::status(std::vector <double> &coordinates) : Rot(coordinates[9], coordinates[10], coordinates[11]), fromSvToA(0,0,0), GOST4401()/*, Damp()*/{

	//скорости
	parametr.vect[0] = coordinates[0];
	parametr.vect[1] = coordinates[1];
	parametr.vect[2] = coordinates[2];
	//кординаты
	parametr.vect[3] = coordinates[3];
	parametr.vect[4] = coordinates[4];
	parametr.vect[5] = coordinates[5];
	//угловые скорости
	parametr.vect[6] = coordinates[6];
	parametr.vect[7] = coordinates[7];
	parametr.vect[8] = coordinates[8];
	//параметры Родрига-Гамильтона
	Vect<4> a = Rot.RG.getRGPar();
	parametr.vect[9] = a.vect[0];
	parametr.vect[10] = a.vect[1];
	parametr.vect[11] = a.vect[2];
	parametr.vect[12] = a.vect[3];

	parametr.vect[13] = coordinates[12];
	parametr.vect[14] = coordinates[13];
	Ve = 0;
	d1 = 0;
	a1 = 0;
	p1 = 0;
	Ve = OMEGA_EARTH * RA_EL*(1 - ALPHA_EL * sin(lat) * sin(lat)) * cos(lat);
	azimut = AZIM;

}

status::~status() {}
void status::addErotation(double lat, double azim) {

	Ve = OMEGA_EARTH * RA_EL*(1 - ALPHA_EL * sin(lat) * sin(lat)) * cos(lat);

	azimut = azim;
	parametr.vect[0] += Ve * sin(azimut);
	parametr.vect[2] += Ve * cos(azimut);


}
//значения производных
Vect<15> status::rightPart() {
	Vect<15> prir;

	return(prir);

}
double status::getH() {
	return h;
}
double status::getV()
{
	return vv;
}
void status::setStageParam(double mDry, double mFuel, double mp, double sm, double sa, double l, double p, double ix, double iy, double iz, double izy) {

	M_DRY = mDry;
	M_FUEL = mFuel;
	M_P = mp;

	S_M = sm;
	S_A = sa;
	L = l;
	P = p;

	I_X = ix;
	I_Y = iy;
	I_Z = iz;
	I_ZY = izy;
}

void status::addV(double Vx, double Vy, double Vz)
{
	parametr.vect[0] += Vx;
	parametr.vect[1] += Vy;
	parametr.vect[2] += Vz;
}

void status::nonIntegr() {

}

//вывод параметров
void status::printParam(std::ofstream &fout) {
	fout << parametr.vect[13] << '\t' << std::scientific;
	fout << parametr.vect[14] << '\t' << std::scientific;
	fout << h << '\t';
	fout << parametr.vect[3] << '\t' << parametr.vect[4] << '\t' << parametr.vect[5] << '\t';
	fout << parametr.vect[0] << '\t' << parametr.vect[1] << '\t' << parametr.vect[2] << '\t';
	fout << sqrt(vFullsq) << '\t';
	fout << parametr.vect[6] << '\t' << parametr.vect[7] << '\t' << parametr.vect[8] << '\t';

	for (int i = 0; i < 3; i++) {
		fout << Rot.Angles.vect[i] * 180 / PI << '\t';
	}

	fout << alpha * 180 / PI << '\t' << betta * 180 / PI << '\t';
	fout << mach << '\t';
	fout << q << '\t';
	for (int i = 0; i < 3; i++) {
		fout << ForcePr.vect[i] << '\t';
	}
	fout << cx << '\t' << cy << '\t' << cz << '\t';
	

	for (int i = 0; i < 3; i++) {
		fout << Torque.vect[i] << '\t';
	}
	for (int i = 0; i < 3; i++) {
		fout << v.vect[i] << '\t';
	}
	fout << mzAl << '\t';
	fout << mzwz << '\t';
	fout << - g * parametr.vect[3] << '\t';
	fout << - g * parametr.vect[4] << '\t';
	fout << - g * parametr.vect[5] << '\t';

	for (int i = 0; i < 3; i++) {
		fout << Fg.vect[i] << '\t';
	}
	fout << Fg.vect[0]+g * parametr.vect[3] << '\t';
	fout << Fg.vect[1]+g * parametr.vect[4] << '\t';
	fout << Fg.vect[2]+g * parametr.vect[5] << '\t';
	fout << myBet << '\t';
	fout << mx << '\t';
	fout << density << '\t' << d1 << '\t';
	fout << pressure << '\t' << p1 << '\t';
	fout << mach << '\t' << a1 << '\t';
	fout << al1 << '\t' << fi1 << '\t';

	fout << Fa.vect[0] << '\t';
	fout << Fa.vect[1] << '\t';
	fout << Fa.vect[2] << '\t';
	Vect<3> FaX = { Fa.vect[0], 0,0 };
	Vect<3> d = FaX * v;
	fout << d.vect[0]* d.vect[0]+ d.vect[1] * d.vect[1] + d.vect[2] * d.vect[2] << '\t';
	double a360 = (alpha < 0) ? 2 * PI + alpha : alpha;
	fout << a360*toDeg << '\t';
	
	fout << mRoll << '\t';
	fout << mPitch << '\t';
	fout << mYaw << '\t';

	fout << '\n';
}
void status::setParam(Vect<15> b) {
	parametr = b;
}

Vect<15> status::getParam() {
	Vect<15> get = parametr;
	return(get);
}

status& status::operator=(const status& right) {
	//проверка на самоприсваивание
	if (this == &right) {
		return *this;
	}
	parametr = right.parametr;
	nonIntegr();

	Fg = right.Fg;
	Torque = right.Torque;

	Rot = right.Rot;

	alpha = right.alpha;
	betta = right.betta;

	g = right.g;
	//ADkoef = right.ADkoef;

	density = right.density;

	m = right.m;//масса
	I_X = right.I_X;
	I_Y = right.I_Y;
	I_Z = right.I_Z;
	I_ZY = right.I_ZY;
	S_M = right.S_M;
	/////TODO: написать модули для определения параметров
	g = right.g;
	/*double Cx = 0.5;
	double Cy = 0.5;
	double Cz = 0.5;*/
	//std::vector <double> ADkoef;
	density = right.density;
	q = right.q;
	mach = right.mach;
	vFullsq = right.vFullsq;
	P = right.P;
	M_P = right.M_P;

	return *this;
}


