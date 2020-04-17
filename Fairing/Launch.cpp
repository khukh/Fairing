#include "pch.h"
#include "Launch.h"
#include "KoefCalc.h"


Launch::Launch(std::vector <double> &coordinates) :status(coordinates)
{
}

Vect<15> Launch::rightPart()
{
	Vect<15> prir;

	//dV/dt
	prir.vect[0] = (Fg.vect[0]) / parametr.vect[14];
	prir.vect[1] = (Fg.vect[1]) / parametr.vect[14];
	prir.vect[2] = (Fg.vect[2]) / parametr.vect[14];
	//dX/dt
	prir.vect[3] = parametr.vect[0];
	prir.vect[4] = parametr.vect[1];
	prir.vect[5] = parametr.vect[2];
	//dW/dt
	prir.vect[6] = Torque.vect[0] / I_X - (I_Z - I_Y) / I_X * parametr.vect[7] * parametr.vect[8];//
	prir.vect[7] = Torque.vect[1] / I_Y - (I_X - I_Z) / I_Y * parametr.vect[6] * parametr.vect[8];
	prir.vect[8] = Torque.vect[2] / I_Z - (I_Y - I_X) / I_Z * parametr.vect[6] * parametr.vect[7];
	//dPRG/dt
	prir.vect[9] = 0;
	prir.vect[10] = 0;
	prir.vect[11] = 0;
	prir.vect[12] = 0;
	//dt/dt
	prir.vect[13] = 1;
	prir.vect[14] = -M_P;
	return prir;
}

void Launch::nonIntegr()
{
	Rot.Angles.vect[0] = pitchProgram();
	//Rot.Angles[0] = 89 * toRad;
	Rot.Angles.vect[1] = 0*toRad;
	Rot.Angles.vect[2] = 0*toRad;
	Rot.fromAnglesToRG();
	
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
		//{ √ринвические координаты }

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



	Ve = OMEGA_EARTH * RA_EL*(1 - ALPHA_EL * sin(LAT) * sin(LAT)) * cos(LAT);

	//std::vector <double> v(3);
	Vect<3> vg = { parametr.vect[0] - Ve * sin(azimut),parametr.vect[1],parametr.vect[2]- Ve * cos(azimut) };
	Rot.fromRGtoMatrixT();
	v = Rot.A*vg;

	vFullsq = v.vect[0] * v.vect[0] + v.vect[1] * v.vect[1] + v.vect[2] * v.vect[2];
//	double r = sqrt(parametr[3] * parametr[3] + parametr[4] * parametr[4] + parametr[5] * parametr[5]);
	vv = sqrt(vFullsq);

	//h = r - RA_EL * (1 - ALPHA_EL * sin(LAT)*sin(LAT));
	density = GOST4401.roFunc(h);
	pressure = GOST4401.pFunc(h);
	double ah = atan2(-v.vect[1], v.vect[0]);
	alpha = ah;
	betta = (vFullsq < 1E-7) ? 0 : asin(v.vect[2]/ sqrt(vFullsq));
	double alphaSpace = sqrt(alpha*alpha + betta * betta);
	mach = vv / GOST4401.aFunc(h);

	q = density * vFullsq / 2;

	

	
	//double height = r - RA_EL * (1 - ALPHA_EL * sin(LAT)*sin(LAT));
	//cx = Cx(mach, alpha);
	//cy = CyAl(mach, alpha, h)*alpha;
		
	//cz = CzBetta(mach, betta, h)*betta;
	double mzwz = MzOmegaZ(mach,alpha);
	//mzAl = MzAlpha(mach, alpha, h);
	//double mzBet = MzAlpha(mach, betta, h);
	//TODO т€га, долгота в радиусе эллипса, пересчитывать азимут, долготу, широту
	Vect<3> ForcePr;
	ForcePr.vect[0] = P - pressure*S_A - cx * density * vFullsq * S_M / 2;
	ForcePr.vect[1] = cy * density * vFullsq * S_M / 2;
	ForcePr.vect[2] = -cz * density * vFullsq * S_M / 2;
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
		Torque.vect[0] = (0 * parametr.vect[6] * L / vv)* density * vFullsq * S_M * L / 2 /*+ Mstab*/;
		//Torque.vect[1] = (0 * parametr.vect[7] * L / vv + mzBet * betta) * density * vFullsq * S_M * L / 2;
		Torque.vect[1] = 0;
		Torque.vect[2] = (mzwz * parametr.vect[8] * L / vv + mzAl * alpha) * density * vFullsq * S_M * L / 2;

	}
	


}

double Launch::pitchProgram()
{
	double pitch = 0;

	double c11 = (39 - 90) * toRad/(67-15);
	double c12 = (27 - 39) * toRad / (118 - 67);
	double c13 = (6 - 27) * toRad / (333 - 118);
	double c14 = (4 - 6) * toRad / (351 - 333);
	double c15 = (-3 - 4) * toRad / (572 - 351);

	double c21 = 1.82756111098252;
	double c22 = 0.955824431533365;
	double c23 = 0.672398241593908;
	double c24 = 0.750491578357562;
	double c25 = 0.263852716330907;



	if ( (parametr.vect[13] >= 0) &&  (parametr.vect[13] < 15) ) { pitch = PI / 2; }
	if ( (parametr.vect[13] >= 15) && (parametr.vect[13] < 67) ) { pitch = c11 * parametr.vect[13] + c21; }
	if ((parametr.vect[13] >= 67) && (parametr.vect[13] < 118) ) { pitch = c12 * parametr.vect[13] + c22; }
	if ((parametr.vect[13] >= 118) && (parametr.vect[13] < 333)) { pitch = c13 * parametr.vect[13] + c23; }														 																				
	if ((parametr.vect[13] >= 333) && (parametr.vect[13] < 351)) { pitch = c14 * parametr.vect[13] + c24; }
	if ((parametr.vect[13] >= 351)&&(parametr.vect[13]<= 572+1)) { pitch = c15 * parametr.vect[13] + c25; }
	return pitch;
}


Launch::~Launch()
{
}
