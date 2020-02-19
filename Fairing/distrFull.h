#pragma once
#include "Drop.h"
class distrFull :
	public Drop
{
public:
	distrFull(std::vector <double> &coordinates, double wx, double wz);
	~distrFull();

	void nonIntegr();

	
	std::random_device rd1;
	std::mt19937 gen1;
	std::normal_distribution<> dis1;
};

