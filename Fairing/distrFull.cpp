#include "pch.h"
#include "distrFull.h"


distrFull::distrFull(std::vector <double> &coordinates, double wx, double wz) :Drop(coordinates), gen1(rd1())
{
	Wx = wx;
	Wz = wz;

	parametr.vect[0] += Wx;
	parametr.vect[2] += Wz;
}


distrFull::~distrFull()
{
}

void distrFull::nonIntegr()
{


	Rot.RG.setRGPar(parametr.vect[9], parametr.vect[10], parametr.vect[11], parametr.vect[12]);
	Rot.RG.norm();
	Rot.fromRGtoAngles();

	double a11 = -cos(LON)*sin(LAT)*cos(AZIM) - sin(LON)*sin(AZIM);
	double a12 = cos(LON)*cos(LAT);
	double a13 = cos(LON)*sin(LAT)*sin(AZIM) - sin(LON)*cos(AZIM);
	double a21 = -sin(LON)*sin(LAT)*cos(AZIM) + cos(LON)*sin(AZIM);
	double a22 = sin(LON)*cos(LAT);
	double a23 = sin(LON)*sin(LAT)*sin(AZIM) + cos(LON)*cos(AZIM);
	double a31 = cos(LAT)*cos(AZIM);
	double a32 = sin(LAT);
	double a33 = -cos(LAT)*sin(AZIM);

	//perevorot: = 1; //0,2,4...тангаж головой вперед
				  //1,3,5...тангаж кормой вперед
		//{ Гринвические координаты }

	double Xgr = parametr.vect[3] * a11 + parametr.vect[4] * a12 + parametr.vect[5] * a13;
	double Ygr = parametr.vect[3] * a21 + parametr.vect[4] * a22 + parametr.vect[5] * a23;
	double Zgr = parametr.vect[3] * a31 + parametr.vect[4] * a32 + parametr.vect[5] * a33;
	double r = sqrt(parametr.vect[3] * parametr.vect[3] + parametr.vect[4] * parametr.vect[4] + parametr.vect[5] * parametr.vect[5]);

	lat = asin(Zgr / r);
	double lat_g = lat * 180 / PI + 0.1921*sin(2 * lat);
	lon = atan2((Ygr / r * cos(lat)), (Xgr / r * cos(lat)));
	double long_C = lon - 2 * PI / 86400 * parametr.vect[13];

	if (long_C > PI) { long_C = -PI + (long_C - PI); }
	if (long_C < -PI) { long_C = PI - abs(long_C + PI); }


	h = r - RA_EL * (1 - ALPHA_EL * sin(lat)*sin(lat));


	//parametr[0] += Wx;
	//parametr[2] += Wz;
	//Ve = 0;
	Ve = OMEGA_EARTH * RA_EL*(1 - ALPHA_EL * sin(LAT) * sin(LAT)) * cos(LAT);
	//std::vector <double> v(3);
	Vect<3> vg = { parametr.vect[0] - Ve * sin(azimut) /*- Wx*/,parametr.vect[1],parametr.vect[2] - Ve * cos(azimut) /*-Wz*/};
	Rot.fromRGtoMatrixT();
	v = Rot.A*vg;

	vFullsq = v.vect[0] * v.vect[0] + v.vect[1] * v.vect[1] + v.vect[2] * v.vect[2];
	vv = sqrt(vFullsq);
//	try {
	//h = 1 / (3-3);
	
	density = GOST4401.roFunc(h);
	d1 = density;
	double dens1 = density*1E6;
//	double qdfg = sqrt(0.05*dens1) / 3;
	std::normal_distribution<> dDensity{ dens1, /*sqrt(*/0.05*dens1/3 };

	density = dDensity(gen1) / 1E6;
	
	pressure = GOST4401.pFunc(h);
	std::normal_distribution<> dPressure{ pressure, /*sqrt(*/0.05*pressure/3 };
	p1 = pressure;
	pressure = dPressure(gen1);
	//pressure *= 1.05;
	double ah = atan2(-v.vect[1], v.vect[0]);
	alpha = ah;
	betta = (vFullsq < 1E-7) ? 0 : atan2(v.vect[2], sqrt(v.vect[0] * v.vect[0] + v.vect[1] * v.vect[1]));
	double alphaSpace = sqrt(alpha*alpha + betta * betta);
	double temperature = pressure / (287.05287*density);
	double aSonic = 20.046796*sqrt(temperature);
	//a1 = vv/ GOST4401.aFunc(h);
	mach = vv / aSonic;

	q = density * vFullsq / 2;

	//ForcePr(3);


	//double height = r - RA_EL * (1 - ALPHA_EL * sin(LAT)*sin(LAT));
	cx = CxPas(mach, alpha, h);
	cy = CyAlPas(mach, alpha, h);
	double signVy = v.vect[1] / abs(v.vect[1]);
	al1 = atan2(-signVy * sqrt(v.vect[1] * v.vect[1] + v.vect[2] * v.vect[2]), v.vect[0]);
	al1 = (al1 < 0) ? 2 * PI + al1 : al1;
	fi1 = atan2(-v.vect[2], -v.vect[1]);
	cz = CzBettaPas(mach, al1, fi1);
	mzwz = MzOmegaZPas(mach, al1, h);
	mzAl = MzAlphaPas(mach, alpha, h);
	myBet = MyBettaPas(mach, al1, fi1);
	mx = MxBettaPas(mach, al1, fi1);

	ForcePr.vect[0] = -cx * q * S_M;
	ForcePr.vect[1] = cy * q * S_M;
	ForcePr.vect[2] = cz * q * S_M;
	Rot.fromRGtoMatrix();
	Fg = Rot.A * ForcePr;

	g = PI0 * parametr.vect[14] / ((r)*(r)*(r));
	Fg.vect[0] -= g * parametr.vect[3];
	Fg.vect[1] -= g * parametr.vect[4];
	Fg.vect[2] -= g * parametr.vect[5];


	if (vFullsq < 1E-5) {
		Torque.vect[0] = 0;
		Torque.vect[1] = 0;
		Torque.vect[2] = 0;

	}
	else {
		double mywy = 0/*MzOmegaZPas(mach, betta, h)*/;
	//	double d = (abs(betta) < 1E-7) ? 0 : betta / abs(betta);
		Torque.vect[0] = (mx)* q * S_M * L /*+ Mstab*/;
		Torque.vect[1] = /*(abs(betta) <PI/2) ?*/ ((mywy * parametr.vect[7] * L / vv + myBet) * q * S_M * L) /*: ((0 * parametr[7] * L / vv - d * mzBet) * density * vFullsq * S_M * L / 2)*/;
		//Torque[1] = 0;
		Torque.vect[2] = /*(abs(alpha)<PI/2)?*/((mzwz * parametr.vect[8] * L / vv + mzAl) * q * S_M * L)/*:((mzwz * parametr[8] * L / vv - s * mzAl) * density * vFullsq * S_M * L / 2)*/;

	}

}
