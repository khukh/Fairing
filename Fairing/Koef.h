#pragma once
class DampingKoef{
public:
	DampingKoef();
	
	~DampingKoef();
	double linInterp(double x1, double x2, double y1, double y2, double x);
	double MrFr(double mach, double alpha);
	double MpFr(double mach, double alpha);
	double MyFr(double mach, double alpha);

	double MrFp(double mach, double alpha);
	double MpFp(double mach, double alpha);
	double MyFp(double mach, double alpha);

	double MrFy(double mach, double alpha);
	double MpFy(double mach, double alpha);
	double MyFy(double mach, double alpha);

	double MrB(double mach, double beta);
	double MpB(double mach, double beta);
	double MyB(double mach, double beta);
};
