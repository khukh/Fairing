#pragma once
#include "status.h"
#include "Koef.h"
class Drop :
	public status
{
public:
	Drop(std::vector <double> &coordinates);
	Vect<15> rightPart();
	
	double Wx, Wz;
	virtual void nonIntegr();
	~Drop();
};

