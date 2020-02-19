#pragma once
#include "status.h"
class Drop :
	public status
{
public:
	Drop(std::vector <double> &coordinates);
	Vect<15> rightPart();
	double mzwz;
	double Wx, Wz;
	virtual void nonIntegr();
	~Drop();
};

