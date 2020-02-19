#pragma once
#include "status.h"
class Launch :
	public status
{
public:
	
	Launch(std::vector<double> &coordinates);
	Vect<15> rightPart();
	void nonIntegr();
	double pitchProgram();
	~Launch();
};

