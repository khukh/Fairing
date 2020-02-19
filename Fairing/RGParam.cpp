#include "pch.h"
#include "RGParam.h"
//#include "matrix.cpp"

RGParam::RGParam(double pitch, double yaw, double roll) {
	//RGPar.resize(4);
	//RGAngle.resize(3);
	/*RGAngle[0] = pitch;
	RGAngle[1] = yaw;
	RGAngle[2] = roll;*/

	RGPar.vect[0] = cos(0.5*yaw)*cos(0.5*pitch)*cos(0.5*roll) - sin(0.5*yaw)*sin(0.5*pitch)*sin(0.5*roll);  //ro
	RGPar.vect[1] = sin(0.5*yaw)*sin(0.5*pitch)*cos(0.5*roll) + cos(0.5*yaw)*cos(0.5*pitch)*sin(0.5*roll);  //l
	RGPar.vect[2] = sin(0.5*yaw)*cos(0.5*pitch)*cos(0.5*roll) + cos(0.5*yaw)*sin(0.5*pitch)*sin(0.5*roll);  //mu
	RGPar.vect[3] = cos(0.5*yaw)*sin(0.5*pitch)*cos(0.5*roll) - sin(0.5*yaw)*cos(0.5*pitch)*sin(0.5*roll);  //nu

	
}


void RGParam::norm() {
	double abs = sqrt(pow(RGPar.vect[0], 2) + pow(RGPar.vect[1], 2) + pow(RGPar.vect[2], 2) + pow(RGPar.vect[3], 2));
	for (int i = 0; i < 4; i++) {
		RGPar.vect[i] /= abs;
	}
}

Vect<4> RGParam::getRGPar() {
	//std::vector <double> rgPar = RGPar;
	//?
	return RGPar;
}

void RGParam::setRGPar(double ro, double lymbda, double mu, double nu) {
	RGPar.vect[0] = ro;
	RGPar.vect[1] = lymbda;
	RGPar.vect[2] = mu;
	RGPar.vect[3] = nu;
}


RGParam::~RGParam() {}
