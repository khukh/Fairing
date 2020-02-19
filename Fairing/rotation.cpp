#include "pch.h"
#include "rotation.h"


rotation::rotation(double pitch, double yaw, double roll) :A(pitch, yaw, roll), RG(pitch, yaw, roll) {
	//Angles.resize(3);
	Angles.vect[0] = pitch;
	Angles.vect[1] = yaw;
	Angles.vect[2] = roll;
}

void rotation::fromAnglesToRG() {
	RG.RGPar.vect[0] = cos(0.5*Angles.vect[1])*cos(0.5*Angles.vect[0])*cos(0.5*Angles.vect[2]) - sin(0.5*Angles.vect[1])*sin(0.5*Angles.vect[0])*sin(0.5*Angles.vect[2]);  //ro
	RG.RGPar.vect[1] = sin(0.5*Angles.vect[1])*sin(0.5*Angles.vect[0])*cos(0.5*Angles.vect[2]) + cos(0.5*Angles.vect[1])*cos(0.5*Angles.vect[0])*sin(0.5*Angles.vect[2]);  //l
	RG.RGPar.vect[2] = sin(0.5*Angles.vect[1])*cos(0.5*Angles.vect[0])*cos(0.5*Angles.vect[2]) + cos(0.5*Angles.vect[1])*sin(0.5*Angles.vect[0])*sin(0.5*Angles.vect[2]);  //mu
	RG.RGPar.vect[3] = cos(0.5*Angles.vect[1])*sin(0.5*Angles.vect[0])*cos(0.5*Angles.vect[2]) - sin(0.5*Angles.vect[1])*cos(0.5*Angles.vect[0])*sin(0.5*Angles.vect[2]);  //nu
}

void rotation::fromRGtoAngles() {
	Angles.vect[0] = asin(2 * (RG.RGPar.vect[0] * RG.RGPar.vect[3] + RG.RGPar.vect[1] * RG.RGPar.vect[2]));
	Angles.vect[2] = atan2(2 * (RG.RGPar.vect[0] * RG.RGPar.vect[1] - RG.RGPar.vect[3] * RG.RGPar.vect[2]), (RG.RGPar.vect[0]* RG.RGPar.vect[0] + RG.RGPar.vect[2]* RG.RGPar.vect[2] - RG.RGPar.vect[3]* RG.RGPar.vect[3] - RG.RGPar.vect[1]* RG.RGPar.vect[1]));
	Angles.vect[1] = atan2(2 * (RG.RGPar.vect[0] * RG.RGPar.vect[2] - RG.RGPar.vect[1] * RG.RGPar.vect[3]), (RG.RGPar.vect[0]* RG.RGPar.vect[0] + RG.RGPar.vect[1]* RG.RGPar.vect[1] - RG.RGPar.vect[3]* RG.RGPar.vect[3] - RG.RGPar.vect[2]* RG.RGPar.vect[2]));
}

void rotation::fromRGtoMatrix() {
	A.matr[0][0] = (RG.RGPar.vect[0]* RG.RGPar.vect[0] + RG.RGPar.vect[1]*RG.RGPar.vect[1] - RG.RGPar.vect[2]* RG.RGPar.vect[2] - RG.RGPar.vect[3]* RG.RGPar.vect[3]);
	A.matr[1][0] = 2*(RG.RGPar.vect[0] * RG.RGPar.vect[3] + RG.RGPar.vect[1] * RG.RGPar.vect[2]);
	A.matr[2][0] = 2 * (-RG.RGPar.vect[0] * RG.RGPar.vect[2] + RG.RGPar.vect[1] * RG.RGPar.vect[3]);

	A.matr[0][1] = 2 * (-RG.RGPar.vect[0] * RG.RGPar.vect[3] + RG.RGPar.vect[1] * RG.RGPar.vect[2]);
	A.matr[1][1] = (RG.RGPar.vect[0] * RG.RGPar.vect[0] + RG.RGPar.vect[2] * RG.RGPar.vect[2] - RG.RGPar.vect[3] * RG.RGPar.vect[3] - RG.RGPar.vect[1] * RG.RGPar.vect[1]);
	A.matr[2][1] = 2 * (RG.RGPar.vect[0] * RG.RGPar.vect[1] + RG.RGPar.vect[3] * RG.RGPar.vect[2]);

	A.matr[0][2] = 2 * (RG.RGPar.vect[0] * RG.RGPar.vect[2] + RG.RGPar.vect[1] * RG.RGPar.vect[3]);
	A.matr[1][2] = 2 * (-RG.RGPar.vect[0] * RG.RGPar.vect[1] + RG.RGPar.vect[2] * RG.RGPar.vect[3]);
	A.matr[2][2] = (RG.RGPar.vect[0] * RG.RGPar.vect[0] + RG.RGPar.vect[3] * RG.RGPar.vect[3] - RG.RGPar.vect[1] * RG.RGPar.vect[1] - RG.RGPar.vect[2] * RG.RGPar.vect[2]);
}

void rotation::fromRGtoMatrixT() {
/*	A.matr[0][0] = (RG.RGPar[0] * RG.RGPar[0] + RG.RGPar[1] * RG.RGPar[1] - RG.RGPar[2] * RG.RGPar[2] - RG.RGPar[3] * RG.RGPar[3]);
	A.matr[0][1] = 2 * (RG.RGPar[0] * RG.RGPar[3] + RG.RGPar[1] * RG.RGPar[2]);
	A.matr[0][2] = 2 * (-RG.RGPar[0] * RG.RGPar[2] + RG.RGPar[1] * RG.RGPar[3]);

	A.matr[1][0] = 2 * (-RG.RGPar[0] * RG.RGPar[3] + RG.RGPar[1] * RG.RGPar[2]);
	A.matr[1][1] = (RG.RGPar[0] * RG.RGPar[0] + RG.RGPar[2] * RG.RGPar[2] - RG.RGPar[3] * RG.RGPar[3] - RG.RGPar[1] * RG.RGPar[1]);
	A.matr[1][2] = 2 * (RG.RGPar[0] * RG.RGPar[1] + RG.RGPar[3] * RG.RGPar[2]);

	A.matr[2][0] = 2 * (RG.RGPar[0] * RG.RGPar[2] + RG.RGPar[1] * RG.RGPar[3]);
	A.matr[2][1] = 2 * (-RG.RGPar[0] * RG.RGPar[1] + RG.RGPar[2] * RG.RGPar[3]);
	A.matr[2][2] = (RG.RGPar[0] * RG.RGPar[0] + RG.RGPar[3] * RG.RGPar[3] - RG.RGPar[1] * RG.RGPar[1] - RG.RGPar[2] * RG.RGPar[2]);*/

	for (int i = 0; i < 3; i++) {
		for (int j = i+1; j < 3; j++) {
			double buf = A.matr[i][j];
			A.matr[i][j] = A.matr[j][i];
			A.matr[j][i] = buf;
		}
		
	}
}


rotation::~rotation() {}
