#pragma once
#include "pch.h"

class rotation {
public:
	matrix A;
	Vect<3> Angles;
	RGParam RG;

	rotation(double pitch, double yaw, double roll);
	void fromAnglesToRG();	
	void fromRGtoAngles();
	void fromRGtoMatrix();
	void fromRGtoMatrixT();
	~rotation();
};

