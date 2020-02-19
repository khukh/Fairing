#pragma once
#include "pch.h"
class RGParam {
public:
	RGParam(double pitch, double yaw, double roll);
	Vect<4> RGPar;
	//std::vector <double> RGAngle;

	void norm();
	Vect<4> getRGPar();
	void setRGPar(double ro, double lymbda, double mu, double nu);
	~RGParam();
};

