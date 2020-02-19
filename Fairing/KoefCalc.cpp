#include "pch.h"
#include "KoefCalc.h"


/////////////////////
double CyAl(double mach, double alpha, double Hg) {
	double Cy_alf = 0;

	Cy_alf = 0.1859 * abs(alpha) * 180 / PI - 0.0223;
	Cy_alf = (abs(alpha) > 1E-8) ? alpha/abs(alpha):0;
	/*if (abs(alpha) == 4.5 * PI / 180) { Cy_alf = 0; }
	if ((abs(alpha) > 0 * PI / 180) and (abs(alpha)) < 4.5*PI / 180) {
		if ((mach >= 0) && (mach <= 0.301)) { Cy_alf = 2.26782; }
		if ((mach >= 0.301) && (mach <= 0.954)) { Cy_alf = 0.525 *mach*mach - 0.624*mach + 2.408; }
		if ((mach >= 0.954) && (mach <= 1.695)) { Cy_alf = -0.718*mach*mach + 2.679*mach + 0.392; }
		if ((mach >= 1.695) && (mach <= 5.613)) { Cy_alf = 0.002 *mach*mach*mach - 0.027*mach*mach + 0.007*mach + 2.948; }
		if ((mach >= 5.613)) { Cy_alf = 2.58773; }
	}
			
	if ((abs(alpha) >= 4.5*PI / 180) && (abs(alpha)) < 7.5*PI / 180) {

		if ((mach >= 0) && (mach <= 0.301)) { Cy_alf = 2.3424; }
		if ((mach >= 0.301) &&(mach <= 0.954)) { Cy_alf = 0.606 *mach*mach - 0.736*mach + 2.509;}
		if ((mach >= 0.954) &&(mach <= 2.033)) { Cy_alf = -0.379*mach*mach + 1.756*mach + 1.041;}
		if ((mach >= 2.033) &&(mach <= 5.613)) { Cy_alf = -0.007*mach*mach + 0.013*mach + 3.013;}
		if (mach >= 5.613) { Cy_alf = 2.84679; }
	}
	if ((abs(alpha) >= 7.5*PI / 180) && (abs(alpha) < 10.5*PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { Cy_alf = 2.5531; }
		if ((mach >= 0.301) && (mach <= 0.954)) { Cy_alf = 0.293 *mach*mach - 0.119*mach + 2.264;}
		if ((mach >= 0.954) && (mach <= 2.033)) { Cy_alf = -0.539*mach*mach + 2.294*mach + 0.759;}
		if ((mach >= 2.033) && (mach <= 4.309)) { Cy_alf = -0.025*mach*mach + 0.290*mach + 2.699;}
		if ((mach >= 4.309) && (mach <= 5.613)) { Cy_alf = -0.075*mach + 3.806; }
		if (mach >= 5.613) { Cy_alf = 3.38318; }
	}
		
	if ((abs(alpha) >= 10.5*PI / 180) && (abs(alpha) <= 12 * PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { Cy_alf = 2.27335; }
		if ((mach >= 0.301) && (mach <= 1.502)) { Cy_alf = 0.787 *mach*mach - 0.604*mach + 2.383; }
		if ((mach >= 1.502) && (mach <= 1.695)) { Cy_alf = -0.515*mach + 4.019; }
		if ((mach >= 1.695) && (mach <= 4.309)) { Cy_alf = 0.076* mach*mach*mach - 0.928*mach*mach + 3.708*mach - 0.868; }
		if ((mach >= 4.309) && (mach <= 5.613)) { Cy_alf = -0.139*mach + 4.588; }
		if (mach >= 5.613) { Cy_alf = 3.80723; }
	}*/
	
	
	if (Hg > 80000) { Cy_alf = 0; }
	return Cy_alf;

}

double CyAlPas(double mach, double alpha, double Hg) {
	double Cy_alf = 0;

	double alf_param = alpha;
	double M = mach;

	if (alf_param > (PI / 2)) { alf_param = abs(alf_param - PI); }
	if (alf_param < (-PI / 2)) { alf_param = abs(alf_param + PI); }

	if (abs(alpha) < 50 * toRad) {
		Cy_alf = 0.0002 * alpha * toDeg * alpha * toDeg + 0.0066 * abs(alpha) * toDeg - 0.0108;
	}
	else {
		Cy_alf = -0.0003 * alpha * toDeg * alpha * toDeg + 0.0562 * abs(alpha) * toDeg - 1.2656;
	}
	if (alpha < 0) {
		Cy_alf *= -1;
	}

	/*if (alf_param == 0 * PI / 180) { Cy_alf = 0; }
	if ((abs(alf_param) > 0 * PI / 180) && (abs(alf_param)) < 4.5*PI / 180) {
		if ((M >= 0) && (M <= 0.301)) { Cy_alf = 0.16592 / 2.5; }
		if ((M >= 0.301) && (M <= 0.954)) { Cy_alf = (-0.061*M + 0.191) / 2.5; }
		if ((M >= 0.954) && (M <= 1.136)) { Cy_alf = (0.6*M - 0.451) / 2.5; }
		if ((M >= 1.136) && (M <= 5.613)) { Cy_alf = (-0.027*M + 0.282) / 2.5; }
		if (M >= 5.613) { Cy_alf = 0.14091 / 2.5; }
	}
		
	if ((abs(alf_param) >= 4.5*PI / 180) && (abs(alf_param) < 8 * PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Cy_alf = 0.36064 / 2.5; }
		if ((M >= 0.301) && (M <= 0.954)) {Cy_alf = (-0.058*M + 0.388) / 2.5;}
		if ((M >= 0.954) && (M <= 1.136)) {Cy_alf = (0.903*M - 0.545) / 2.5;}
		if ((M >= 1.136) && (M <= 5.613)) {Cy_alf = (-0.056*M + 0.586) / 2.5; }
		if (M >= 5.613) { Cy_alf = 0.27924 / 2.5; }
	}
		
	
	if ((abs(alf_param) >= 8 * PI / 180) && (abs(alf_param) < 12.5*PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Cy_alf = 0.62224 / 2.5; }
		if ((M >= 0.301) && (M <= 0.954)) { Cy_alf = (0.015*M + 0.618) / 2.5; }
		if ((M >= 0.954) && (M <= 1.136)) { Cy_alf = (1.157*M - 0.473) / 2.5; }
		if ((M >= 1.136) && (M <= 5.613)) { Cy_alf = (0.012*M*M*M - 0.121*M*M + 0.302*M - 0.663) / 2.5; }
		if (M >= 5.613) { Cy_alf = 0.72451 / 2.5; }
	}
		
	if ((abs(alf_param) >= 12.5*PI / 180) && (abs(alf_param) <= 17.5*PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Cy_alf = 0.98213 / 2.5; }
		if ((M >= 0.301) && (M <= 0.954)) { Cy_alf = (0.164*M + 0.933) / 2.5; }
		if ((M >= 0.954) && (M <= 1.316)) { Cy_alf = (-2.552*M*M + 7.025*M - 3.291) / 2.5; }
		if ((M >= 1.316) && (M <= 5.613)) { Cy_alf = (0.002*M*M - 0.154*M + 1.763) / 2.5; }
		if (M >= 5.613) { Cy_alf = 0.99124 / 2.5; }
	}
			
	if ((abs(alf_param) >= 17.5*PI / 180) && (abs(alf_param) <= 25 * PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Cy_alf = 1.39881 / 2.5; }
		if ((M >= 0.301) && (M <= 0.954)) { Cy_alf = (1.154*M*M - 1.103*M + 1.626) / 2.5; }
		if ((M >= 0.954) && (M <= 1.316)) { Cy_alf = (-6.024*M*M + 15.74*M - 7.917) / 2.5; }
		if ((M >= 1.316) && (M <= 5.613)) { Cy_alf = (0.042*M*M - 0.51*M + 3.042) / 2.5; }
		if (M >= 5.613) { Cy_alf = 1.54017 / 2.5; }
	}
			
	if ((abs(alf_param) >= 25 * PI / 180) && (abs(alf_param) <= 37.5*PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Cy_alf = 2.33457 / 2.5; }
		if ((M >= 0.301) && (M <= 1.136)) { Cy_alf = (1.88*M*M - 0.802*M + 2.406) / 2.5; }
		if ((M >= 1.136) && (M <= 1.695)) { Cy_alf = (-1.778*M*M + 4.936*M + 0.62) / 2.5; }
		if ((M >= 1.695) && (M <= 5.613)) { Cy_alf = (0.08*M*M - 0.855*M + 5.098) / 2.5; }
		if (M >= 5.613) { Cy_alf = 2.83184 / 2.5; }
	}

	if ((abs(alf_param) >= 37.5*PI / 180) && (abs(alf_param) <= 52.5*PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Cy_alf = 3.79029 / 2.5; }
		if ((M >= 0.301) && (M <= 0.702)) { Cy_alf = (-0.05*M + 3.805) / 2.5; }
		if ((M >= 0.702) && (M <= 1.316)) { Cy_alf = (-10.29*M*M + 24.8*M - 8.572) / 2.5; }
		if ((M >= 1.316) && (M <= 5.613)) { Cy_alf = (-0.02*M*M*M + 0.317*M*M - 1.663*M + 7.939) / 2.5; }
		if (M >= 5.613) { Cy_alf = 4.96039 / 2.5; }
	}
		

	if ((abs(alf_param) >= 52.5*PI / 180) && (abs(alf_param) <= 70 * PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Cy_alf = (3.35042) / 2.5; }
		if ((M >= 0.301) && (M <= 0.954)) { Cy_alf = (-0.166*M*M + 2.101*M + 2.733) / 2.5; }
		if ((M >= 0.954) && (M <= 1.136)) { Cy_alf = (19.36*M - 13.88) / 2.5; }
		if ((M >= 1.136) && (M <= 1.502)) { Cy_alf = (-5.248*M*M + 13.12*M - 0.029) / 2.5; }
		if ((M >= 1.502) && (M <= 5.613)) { Cy_alf = (-0.034*M*M*M + 0.447*M*M - 1.995*M + 9.937) / 2.5; }
		if (M >= 5.613) { Cy_alf = 6.80484 / 2.5; }
	}
		
	if ((abs(alf_param) >= 70 * PI / 180) && (abs(alf_param) <= 90 * PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Cy_alf = 3.34229 / 2.5; }
		if ((M >= 0.301) && (M <= 0.702)) { Cy_alf = (2.442*M + 2.607) / 2.5; }
		if ((M >= 0.702) && (M <= 1.136)) { Cy_alf = (-8.778*M*M + 26.08*M - 9.663) / 2.5; }
		if ((M >= 1.136) && (M <= 3.366)) { Cy_alf = (-0.123*M*M*M + 0.891*M*M - 2.513*M + 10.56) / 2.5; }
		if ((M >= 3.366) && (M <= 3.995)) { Cy_alf = (-0.186*M*M + 1.049*M + 6.091) / 2.5; }
		if ((M >= 3.995) && (M <= 5.613)) { Cy_alf = (0.065*M + 7.144) / 2.5; }
		if (M >= 5.613) { Cy_alf = 7.48732 / 2.5; }
	}*/

	if (Hg > 80000) { Cy_alf = 0; }
	return Cy_alf;
}



double CzBetta(double mach, double betta, double Hg) {
	return -CyAl(mach, betta, Hg);
}

double CzBettaPas(double mach, double betta, double Hg) {
	return -CyAlPas(mach, betta, Hg);
}
///////////////////////

double Cx(double mach, double alphaSpace) {
	double fCx = 0;
	/*if (mach == 0) { fCx = 0.187; }
	if ((mach > 0) && (mach <= 0.301)) { fCx = -5 * 0.00001*mach + 0.187; }
	if ((mach >= 0.301) && (mach <= 0.702)) { fCx = 0.022*mach + 0.180; }
	if ((mach >= 0.702) && (mach <= 0.850)) { fCx = 0.149*mach + 0.091; }
	if ((mach >= 0.850) && (mach <= 1.136)) { fCx = 1.747*mach*mach - 2.602*mach + 1.167; }
	if ((mach >= 1.136) && (mach <= 5.613)) { fCx = -0.005*mach*mach*mach + 0.078*mach*mach - 0.372*mach + 0.798; }
	if ((mach >= 5.613)) { fCx = 0.147; }*/
	if ((mach > -1E-6)&&(mach < 0.6)) { fCx = 0.0248 * mach + 0.16; }
	else if (mach < 1.05) { fCx = 1.4757 * mach * mach - 1.9241 * mach + 0.7984; }
	else if (mach < 1.6) { fCx = 1.4976*mach*mach*mach -6.5187*mach*mach+9.2434*mach-3.8497; }
	else { fCx = -0.0189 * mach*mach*mach + 0.1685*mach*mach - 0.5222*mach + 0.8602; }
	
	return fCx;
}

double CxPas(double mach, double alpha, double Hg) {


	double fCx = 0;
	/*if (mach == 0) { fCx = 0; }
	if ((mach > 0) and (mach <= 0.301)) { fCx = -5 * 0.00001*mach + 0.187; }
	if ((mach >= 0.301) && (mach <= 0.702)) { fCx = 0.022*mach + 0.180; }
	if ((mach >= 0.702) && (mach <= 0.85)) { fCx = 0.149*mach + 0.091; }
	if ((mach >= 0.85) && (mach <= 1.136)) { fCx = 1.747*mach*mach - 2.602*mach + 1.167; }
	if ((mach >= 1.136) && (mach <= 5.613)) { fCx = -0.005*mach*mach*mach + 0.078*mach*mach - 0.372*mach + 0.798; }
	if (mach >= 5.613) { fCx = 0.147; }*/

	//if (Hg > 80000) { fCx = 0; }


	if (alpha > (PI / 2)) { alpha = abs(alpha - PI); }
	if (alpha < (-PI / 2)) { alpha = abs(alpha + PI); }

	if (abs(alpha) < 62.5*toRad) {
		fCx = -0.0003*alpha*alpha*toDeg*toDeg + 0.0068*abs(alpha)*toDeg + 1.5777;
	}
	else {
		fCx = 0.0006*alpha*alpha*toDeg*toDeg - 0.1242 * abs(alpha)*toDeg + 6.0343;
	}

	/*if ((abs(alpha) >= 0 * PI / 180) && (abs(alpha) < 30 * PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { fCx = (0.137); }
		if ((mach >= 0.301) && (mach <= 0.702)) { fCx = (-0.0063*mach + 0.1389); }
		if ((mach >= 0.702) && (mach <= 0.954)) { fCx = (0.1995*mach - 0.0055); }
		if ((mach >= 0.954) && (mach <= 1.136)) { fCx = (1.1201*mach - 0.8839); }
		if ((mach >= 1.136) && (mach <= 2.372)) { fCx = (-0.0475*mach*mach + 0.1917*mach + 0.2214); }
		if ((mach >= 2.372) && (mach <= 5.613)) { fCx = (-0.0047*mach + 0.416); }
		if (mach >= 5.613) { fCx = 0.37; }
		if (Hg > 80000) { fCx = 0; }

	}

	if ((abs(alpha) >= 30 * PI / 180) && (abs(alpha) < 50 * PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { fCx = (0.291) / 4; }
		if ((mach >= 0.301) && (mach <= 0.954)) { fCx = (0.6693*mach*mach - 0.6844*mach + 0.4366) / 4; }
		if ((mach >= 0.954) && (mach <= 1.136)) { fCx = (2.3807*mach - 1.8787) / 4; }
		if ((mach >= 1.136) && (mach <= 1.502)) { fCx = (1.5525*mach*mach - 3.9963*mach + 3.3619) / 4; }
		if ((mach >= 1.502) && (mach <= 3.366)) { fCx = (0.0662*mach*mach*mach - 0.4865*mach*mach + 1.1339*mach + 0.0331) / 4; }
		if ((mach >= 3.366) && (mach <= 5.613)) { fCx = (-0.012*mach + 0.8946) / 4; }
		if (mach >= 5.613) { fCx = 0.831 / 4; }
	}

	if ((abs(alpha) >= 50 * PI / 180) && (abs(alpha) < 80 * PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { fCx = (0.291) / 3; }
		if ((mach >= 0.301) && (mach <= 0.954)) { fCx = (0.6693*mach*mach - 0.6844*mach + 0.4366) / 3; }
		if ((mach >= 0.954) && (mach <= 1.136)) { fCx = (2.3807*mach - 1.8787) / 3; }
		if ((mach >= 1.136) && (mach <= 1.502)) { fCx = (1.5525*mach*mach - 3.9963*mach + 3.3619) / 3; }
		if ((mach >= 1.502) && (mach <= 3.366)) { fCx = (0.0662*mach*mach*mach - 0.4865*mach*mach + 1.1339*mach + 0.0331) / 3; }
		if ((mach >= 3.366) && (mach <= 5.613)) { fCx = (-0.012*mach + 0.8946) / 3; }
		if (mach >= 5.613) { fCx = 0.831 / 3;	}
	}

	if (abs(alpha) >= 80 * PI / 180) {
		if ((mach >= 0) && (mach <= 0.301)) { fCx = (0.291) / 2; }
		if ((mach >= 0.301) && (mach <= 0.954)) { fCx = (0.6693*mach*mach - 0.6844*mach + 0.4366) / 2; }
		if ((mach >= 0.954) && (mach <= 1.136)) { fCx = (2.3807*mach - 1.8787) / 2; }
		if ((mach >= 1.136) && (mach <= 1.502)) { fCx = (1.5525*mach*mach - 3.9963*mach + 3.3619) / 2; }
if ((mach >= 1.502) && (mach <= 3.366)) { fCx = (0.0662*mach*mach*mach - 0.4865*mach*mach + 1.1339*mach + 0.0331) / 2; }
if ((mach >= 3.366) && (mach <= 5.613)) { fCx = (-0.012*mach + 0.8946) / 2; }
if (mach >= 5.613) { fCx = 0.831 / 2; }
	}

	if ((abs(alpha) >= 0 * PI / 180) && (abs(alpha) < 15 * PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { fCx = (0.137); }
		if ((mach >= 0.301) && (mach <= 0.702)) { fCx = (-0.0063*mach + 0.1389); }
		if ((mach >= 0.702) && (mach <= 0.954)) { fCx = (0.1995*mach - 0.0055); }
		if ((mach >= 0.954) && (mach <= 1.136)) { fCx = (1.1201*mach - 0.8839); }
		if ((mach >= 1.136) && (mach <= 2.372)) { fCx = (-0.0475*mach*mach + 0.1917*mach + 0.2214); }
		if ((mach >= 2.372) && (mach <= 5.613)) { fCx = (-0.0047*mach + 0.416); }
		if (mach >= 5.613) { fCx = 0.37; }
	}

	if ((abs(alpha) >= 15 * PI / 180) && (abs(alpha) < 37.5*PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { fCx = (0.144) / 2.5; }
		if ((mach >= 0.301) && (mach <= 0.954)) { fCx = (0.7126*mach*mach - 0.438*mach + 0.2113) / 2.5; }
		if ((mach >= 0.954) && (mach <= 1.136)) { fCx = (2.3987*mach - 1.8468) / 2.5; }
		if ((mach >= 1.136) && (mach <= 1.695)) { fCx = (-1.1626*mach*mach*mach + 4.771*mach*mach - 6.3684*mach + 3.6601) / 2.5; }
		if ((mach >= 1.695) && (mach <= 5.613)) { fCx = (-0.0028*mach*mach*mach + 0.0393*mach*mach - 0.2033*mach + 1.1585) / 2.5; }
		if (mach >= 5.613) { fCx = 0.769 / 2.5; }
	}
	if ((abs(alpha) >= 37.5*PI / 180) && (abs(alpha) < 52.5*PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { fCx = 0.192 / 2.5; }
		if ((mach >= 0.301) and (mach <= 0.954)) { fCx = (1.0722*mach*mach - 1.3759*mach + 0.5083) / 2.5; }
		if ((mach >= 0.954) and (mach <= 1.136)) { fCx = (3.4154*mach - 3.0874) / 2.5; }
		if ((mach >= 1.136) and (mach <= 1.695)) { fCx = (-1.4862*mach*mach*mach + 6.1585*mach*mach - 8.3367*mach + 4.494) / 2.5; }
		if ((mach >= 1.695) and (mach <= 5.613)) { fCx = (0.0102*mach*mach - 0.0991*mach + 0.9584) / 2.5; }
		if (mach >= 5.613) { fCx = 0.725 / 2.5; }
	}
	if ((abs(alpha) >= 52.5*PI / 180) && (abs(alpha) <= 70 * PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { fCx = 0.281 / 2.5; }
		if ((mach >= 0.301) && (mach <= 0.954)) { fCx = (0.2438*mach*mach - 0.5964*mach + 0.4385) / 2.5; }
		if ((mach >= 0.954) && (mach <= 1.502)) { fCx = (5.5606*mach*mach*mach - 22.752*mach*mach + 31.194*mach - 13.789) / 2.5; }
		if ((mach >= 1.502) && (mach <= 2.033)) { fCx = (-0.1695*mach*mach + 0.6908*mach - 0.0779) / 2.5; }
		if ((mach >= 2.033) && (mach <= 5.613)) { fCx = (0.0058*mach*mach - 0.0602*mach + 0.7301) / 2.5; }
		if (mach >= 5.613) { fCx = 0.578 / 2.5; }
	}
	if ((abs(alpha) >= 70 * PI / 180) && (abs(alpha) <= 90 * PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { fCx = 0.18 / 2.5; }
		if ((mach >= 0.301) && (mach <= 0.702)) { fCx = (-0.275*mach + 0.262) / 2.5; }
		if ((mach >= 0.702) && (mach <= 0.954)) { fCx = (1.452*mach - 0.95) / 2.5; }
		if ((mach >= 0.954) && (mach <= 1.136)) { fCx = (-0.553*mach + 0.963) / 2.5; }
		if ((mach >= 1.136) && (mach <= 1.502)) { fCx = (-0.658*mach*mach + 2.105*mach - 1.206) / 2.5; }
		if ((mach >= 1.502) && (mach <= 2.372)) { fCx = (-0.121*mach*mach + 0.486*mach + 0.012) / 2.5; }
		if ((mach >= 2.372) && (mach <= 4.309)) { fCx = (-0.003*mach*mach + 0.007*mach + 0.485) / 2.5; }
		if ((mach >= 4.309) && (mach <= 5.613)) { fCx = (-0.0007*mach + 0.4621) / 2.5; }
		if (mach >= 5.613) { fCx = 0.4581709 / 2.5; }
	}*/

	if (Hg > 80000) { fCx = 0; }
	return fCx;

}




////////////////////////
double MxOmegaX(double mach, double alpha) {
	return -0.005*0.6786;

}

/////////
double MzOmegaZ(double mach, double alpha) {
	double m = 0;
	return m;

}

double MzOmegaZPas(double mach, double alpha, double Hg) {
	double Mz_wz = 0;
	double M = mach;
	if ((M >= 0) && (M <= 1.136)) { Mz_wz = -1.2174; }
	if ((M >= 1.136) && (M <= 5.613)) { Mz_wz = 0.0142*M*M*M - 0.19*M*M + 0.909*M - 1.9966; }
	if (M >= 5.613) { Mz_wz = -0.382; }
	if (Hg > 80000) { Mz_wz = 0; }
	return Mz_wz;
}

double MyOmegaYPas(double mach, double alpha, double Hg) {
	double Mz_wz = 0;
	double M = mach;
	if ((M >= 0) && (M <= 1.136)) { Mz_wz = -1.2174; }
	if ((M >= 1.136) && (M <= 5.613)) { Mz_wz = 0.0142*M*M*M - 0.19*M*M + 0.909*M - 1.9966; }
	if (M >= 5.613) { Mz_wz = -0.382; }
	if (Hg > 80000) { Mz_wz = 0; }
	return Mz_wz;
}

double MzAlphaPas(double mach, double alpha, double Hg) {
	double Mz_alf = 0;
	double M = mach;
	double alf_param = alpha;
	double xd = 0;
	if (alpha > (PI / 2)) { alpha = abs(alpha - PI); }
	if (alpha < (-PI / 2)) { alpha = abs(alpha + PI); }

	if (abs(alpha) < 6.25*toRad) {
		xd = -0.0008*alpha*alpha*toDeg*toDeg - 0.0011*abs(alpha)*toDeg + 0.2987;
	}
	else {
		if (abs(alpha) < 25 * toRad) {
			xd = 0.0005*alpha*alpha*toDeg*toDeg - 0.0255 * abs(alpha)*toDeg - 0.2033;
		}
		else {
			xd = -0.5;
		}
		
	}

	Mz_alf = xd * CyAlPas(mach, alpha, Hg);

	/*if (alf_param == 0 * PI / 180) { Mz_alf = 0; }
	if ((abs(alf_param) > 0 * PI / 180) && (abs(alf_param) < 4.5*PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Mz_alf = 0.07043 / 2.5 / 1; }
		if ((M >= 0.301) && (M <= 1.316)) { Mz_alf = (-0.086*M*M*M + 0.218*M*M - 0.147*M + 0.097) / 2.5 / 1; }
		if ((M >= 1.316) && (M <= 5.613)) { Mz_alf = (-0.001*M*M*M + 0.021*M*M - 0.101*M + 0.187) / 2.5 / 1; }
		if (M >= 5.613) { Mz_alf = 0.01879 / 2.5 / 1; }
	}
		
	if ((abs(alf_param) >= 4.5*PI / 180) && (abs(alf_param) < 8 * PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Mz_alf = 0.13587 / 2.5 / 1; }
		if ((M >= 0.301) && (M <= 0.954)) { Mz_alf = (0.053*M*M - 0.052*M + 0.146) / 2.5 / 1; }
		if ((M >= 0.954) && (M <= 1.502)) { Mz_alf = (-0.27*M*M + 0.673*M - 0.25) / 2.5 / 1; }
		if ((M >= 1.502) && (M <= 5.613)) { Mz_alf = (-0.002*M*M*M + 0.037*M*M - 0.188*M + 0.357) / 2.5 / 1; }
		if (M >= 5.613) { Mz_alf = 0.03746 / 2.5 / 1; }
	}
		
	if ((abs(alf_param) >= 8 * PI / 180) && (abs(alf_param) < 12.5*PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Mz_alf = 0.2235 / 2.5 / 1; }
		if ((M >= 0.301) && (M <= 0.954)) { Mz_alf = (-0.003*M + 0.224) / 2.5 / 1; }
		if ((M >= 0.954) && (M <= 1.136)) { Mz_alf = (0.274*M - 0.041) / 2.5 / 1; }
		if ((M >= 1.136) && (M <= 2.711)) { Mz_alf = (-0.015*M*M - 0.059*M + 0.365) / 2.5 / 1; }
		if ((M >= 2.711) && (M <= 4.309)) { Mz_alf = (-0.009*M*M + 0.036*M + 0.063) / 2.5 / 1; }
		if ((M >= 4.309) && (M <= 5.613)) { Mz_alf = (-0.003*M + 0.06) / 2.5 / 1; }
		if (M >= 5.613) { Mz_alf = 0.03818 / 2.5 / 1; }
	}
		
	if ((abs(alf_param) >= 12.5*PI / 180) && (abs(alf_param) <= 17.5*PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Mz_alf = 0.31623 / 2.5 / 1; }
		if ((M >= 0.301) and (M <= 0.954)) { Mz_alf = (0.097*M*M - 0.074*M + 0.329) / 2.5 / 1; }
		if ((M >= 0.954) and (M <= 1.136)) { Mz_alf = (0.314*M + 0.047) / 2.5 / 1; }
		if ((M >= 1.136) and (M <= 3.682)) { Mz_alf = (0.034*M*M - 0.291*M + 0.693) / 2.5 / 1; }
		if ((M >= 3.682) and (M <= 5.613)) { Mz_alf = (0.01*M*M - 0.105*M + 0.336) / 2.5 / 1; }
		if (M >= 5.613) { Mz_alf = 0.07561 / 2.5 / 1; }
	}
		
	if ((abs(alf_param) >= 17.5*PI / 180) && (abs(alf_param) <= 25 * PI / 180)) {
		if ((M >= 0) and (M <= 0.301)) { Mz_alf = 0.41185 / 2.5 / 1; }
		if ((M >= 0.301) and (M <= 0.702)) { Mz_alf = (-0.112*M + 0.445) / 2.5 / 1; }
		if ((M >= 0.702) and (M <= 1.316)) { Mz_alf = (-2.855*M*M*M + 8.081*M*M - 7.186*M + 2.416) / 2.5 / 1; }
		if ((M >= 1.316) and (M <= 3.366)) { Mz_alf = (0.064*M*M - 0.44*M + 0.914) / 2.5 / 1; }
		if ((M >= 3.366) and (M <= 5.613)) { Mz_alf = (-0.009*M*M + 0.06*M + 0.065) / 2.5 / 1; }
		if (M >= 5.613) { Mz_alf = 0.10175 / 2.5 / 1; }
	}
		
	if ((abs(alf_param) >= 25 * PI / 180) && (abs(alf_param) <= 37.5*PI / 180)) {
		if ((M >= 0) and (M <= 0.301)) { Mz_alf = 0.5233 / 2.5 / 1; }
		if ((M >= 0.301) and (M <= 0.954)) { Mz_alf = (0.316*M*M - 0.302*M + 0.585) / 2.5 / 1; }
		if ((M >= 0.954) and (M <= 1.316)) { Mz_alf = (-1.711*M*M + 3.824*M - 1.505) / 2.5 / 1; }
		if ((M >= 1.316) and (M <= 2.711)) { Mz_alf = (0.109*M*M - 0.626*M + 1.197) / 2.5 / 1; }
		if ((M >= 2.711) and (M <= 5.613)) { Mz_alf = (0.004*M*M - 0.063*M + 0.438) / 2.5 / 1; }
		if (M >= 5.613) { Mz_alf = 0.2164 / 2.5 / 1; }
	}
		
	if ((abs(alf_param) >= 37.5*PI / 180) && (abs(alf_param) <= 52.5*PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Mz_alf = 070114 / 2.5 / 1; }
		if ((M >= 0.301) && (M <= 0.954)) { Mz_alf = (0.024*M + 0.693) / 2.5 / 1; }
		if ((M >= 0.954) && (M <= 1.136)) { Mz_alf = (0.224*M + 0.502) / 2.5 / 1; }
		if ((M >= 1.136) && (M <= 3.044)) { Mz_alf = (-0.085*M*M*M + 0.674*M*M - 1.81*M + 2.066) / 2.5 / 1; }
		if ((M >= 3.044) && (M <= 5.613)) { Mz_alf = (0.005*M*M - 0.064*M + 0.553) / 2.5 / 1; }
		if (M >= 5.613) { Mz_alf = 0.36198 / 2.5 / 1; }
	}
		
	if ((abs(alf_param) >= 52.5*PI / 180) && (abs(alf_param) <= 70 * PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Mz_alf = 0.4066 / 2.5 / 1; }
		if ((M >= 0.301) && (M <= 0.954)) { Mz_alf = (-0.043*M + 0.419) / 2.5 / 1; }
		if ((M >= 0.954) && (M <= 1.136)) { Mz_alf = (1.891*M - 1.425) / 2.5 / 1; }
		if ((M >= 1.136) && (M <= 1.316)) { Mz_alf = (-0.813*M + 1.645) / 2.5 / 1; }
		if ((M >= 1.316) && (M <= 1.502)) { Mz_alf = (-0.139*M + 0.758) / 2.5 / 1; }
		if ((M >= 1.502) && (M <= 2.711)) { Mz_alf = (0.12*M*M - 0.641*M + 1.237) / 2.5 / 1; }
		if ((M >= 2.711) && (M <= 3.366)) { Mz_alf = (-0.034*M + 0.476) / 2.5 / 1; }
		if ((M >= 3.366) && (M <= 5.613)) { Mz_alf = (0.006*M*M - 0.062*M + 0.497) / 2.5 / 1; }
		if (M >= 5.613) { Mz_alf = 0.34737 / 2.5 / 1; }
	}

	if ((abs(alf_param) >= 70 * PI / 180) && (abs(alf_param) <= 90 * PI / 180)) {
		if ((M >= 0) && (M <= 0.301)) { Mz_alf = 0.29464 / 2.5 / 1; }
		if ((M >= 0.301) && (M <= 0.702)) { Mz_alf = (-0.045*M + 0.308) / 2.5 / 1; }
		if ((M >= 0.702) && (M <= 0.954)) { Mz_alf = (2.034*M - 1.151) / 2.5 / 1; }
		if ((M >= 0.954) && (M <= 1.136)) { Mz_alf = (-1.317*M + 2.046) / 2.5 / 1; }
		if ((M >= 1.136) && (M <= 1.502)) { Mz_alf = (2.115*M*M - 5.703*M + 4.299) / 2.5 / 1; }
		if ((M >= 1.502) && (M <= 2.711)) { Mz_alf = (0.076*M*M - 0.46*M + 1.026) / 2.5 / 1; }
		if ((M >= 2.711) && (M <= 4.309)) { Mz_alf = (-0.026*M*M*M + 0.225*M*M - 0.599*M + 0.833) / 2.5 / 1; }
		if ((M >= 4.309) && (M <= 5.613)) { Mz_alf = (0.003*M + 0.309) / 2.5 / 1; }
		if (M >= 5.613) { Mz_alf = 0.32634 / 2.5 / 1; }
	}*/
		
	if (Hg > 80000) { Mz_alf = 0; }

	return Mz_alf;
}

double MzAlpha(double mach, double alpha, double Hg) {
	double Mz_alf = 0;
	/*if (abs(alpha) == 0 * PI / 180) { Mz_alf = 0; }
	if ((abs(alpha) > 0 * PI / 180) && (abs(alpha) < 4.5*PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { Mz_alf = 0.332204; }
		if ((mach >= 0.301) && (mach <= 0.954)) { Mz_alf = 0.041*mach*mach - 0.112*mach + 0.362; }
		if ((mach >= 0.954) && (mach <= 1.316)) { Mz_alf = 0.880*mach*mach - 2.004*mach + 1.403; }
		if ((mach >= 1.316) && (mach <= 2.711)) { Mz_alf = 0.059*mach*mach*mach - 0.506*mach*mach + 1.313*mach - 0.697; }
		if ((mach >= 2.711) && (mach <= 4.309)) { Mz_alf = 0.0267*mach*mach*mach - 0.2465*mach*mach + 0.7117*mach - 0.315; }
		if ((mach >= 4.309) && (mach <= 5.613)) { Mz_alf = -0.0071*mach + 0.3422; }
		if (mach >= 5.613) { Mz_alf = 0.302512; }
	}
		
	
	if ((abs(alpha) >= 4.5*PI / 180) && (abs(alpha) < 7.5*PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { Mz_alf = 0.38242; }
		if ((mach >= 0.301) & (mach <= 0.954)) { Mz_alf = 0.114*mach*mach - 0.201*mach + 0.432; }
		if ((mach >= 0.954) & (mach <= 1.316)) { Mz_alf = -1.615*mach*mach + 3.681*mach - 1.696; }
		if ((mach >= 1.316) & (mach <= 1.695)) { Mz_alf = 0.262*mach + 0.003; }
		if ((mach >= 1.695) & (mach <= 2.372)) { Mz_alf = -0.289*mach*mach + 1.154*mach - 0.678; }
		if ((mach >= 2.372) & (mach <= 4.309)) { Mz_alf = -0.0207*mach*mach*mach + 0.2356*mach*mach - 0.8821*mach + 1.4773; }
		if ((mach >= 4.309) & (mach <= 5.613)) { Mz_alf = -0.0106*mach + 0.4419; }
		if (mach >= 5.613) { Mz_alf = 0.382609; }
	}

	if ((abs(alpha) >= 7.5*PI / 180) && (abs(alpha) < 10.5*PI / 180)) {
		if ((mach >= 0) && (mach <= 0.301)) { Mz_alf = 0.345199; }
		if ((mach >= 0.301) && (mach <= 0.954)) { Mz_alf = 0.029*mach + 0.337; }
		if ((mach >= 0.954) && (mach <= 1.502)) { Mz_alf = 3.044*mach*mach*mach - 11.28*mach*mach + 13.9*mach - 5.274; }
		if ((mach >= 1.502) && (mach <= 2.711)) { Mz_alf = 0.233*mach*mach*mach - 1.575*mach*mach + 3.44*mach - 1.934; }
		if ((mach >= 2.711) && (mach <= 4.309)) { Mz_alf = 0.0451*mach*mach*mach - 0.4791*mach*mach + 1.7306*mach - 1.5945; }
		if ((mach >= 4.309) && (mach <= 5.613)) { Mz_alf = 0.0105*mach + 0.5313; }
		if (mach >= 5.613) { Mz_alf = 0.590037; }
	
	}
		

	if ((abs(alpha) >= 10.5*PI / 180) && (abs(alpha)) <= 12 * PI / 180) {
		if ((mach >= 0) && (mach <= 0.301)) { Mz_alf = 0.355042; }
		if ((mach >= 0.301) && (mach <= 1.316)) { Mz_alf = 0.155*mach*mach*mach - 0.241*mach*mach + 0.154*mach + 0.326; }
		if ((mach >= 1.316) && (mach <= 1.695)) { Mz_alf = -2.539*mach*mach + 7.669*mach - 5.23; }
		if ((mach >= 1.695) && (mach <= 4.309)) { Mz_alf = 0.0302*mach*mach*mach*mach - 0.3519*mach*mach*mach + 1.4414*mach*mach - 2.3267*mach + 1.7414; }
		if ((mach >= 4.309) && (mach <= 5.613)) { Mz_alf = -0.003*mach + 0.747; }
		if (mach >= 5.613) { Mz_alf = 0.728181; }
	}
		
*/

	Mz_alf = -0.0974 * alpha * 180 / PI + 0.0178 + CyAl(mach, alpha, Hg)*35.04/58;
	if (Hg > 80000) { Mz_alf = 0; }

	return Mz_alf;

}


double MyOmegaY(double mach, double betta) {
	return 0;

}



double MyBettaPas(double mach, double betta, double Hg) {
	return MzAlphaPas(mach, betta, Hg);
}



/*Function Cz_bet(M, bet_param, Hg:real) :real;
begin
Cz_bet : = 0;
if (bet_param = 0 * pi / 180) then Cz_bet : = 0;
if (abs(bet_param) > 0 * pi / 180) and (abs(bet_param) < 4.5*pi / 180) then begin
	if (M >= 0) and (M <= 0.301) then Cz_bet : = -0.16592 / 2.5;
if (M >= 0.301) and (M <= 0.954) then Cz_bet : = -(-0.061*M + 0.191) / 2.5;
if (M >= 0.954) and (M <= 1.136) then Cz_bet : = -(0.6*M - 0.451) / 2.5;
if (M >= 1.136) and (M <= 5.613) then Cz_bet : = -(-0.027*M + 0.282) / 2.5;
if (M >= 5.613) then Cz_bet : = -0.14091 / 2.5;
end;
if (abs(bet_param) >= 4.5*pi / 180) and (abs(bet_param) < 8 * pi / 180) then begin
	if (M >= 0) and (M <= 0.301) then Cz_bet : = -0.36064 / 2.5;
if (M >= 0.301) and (M <= 0.954) then Cz_bet : = -(-0.058*M + 0.388) / 2.5;
if (M >= 0.954) and (M <= 1.136) then Cz_bet : = -(0.903*M - 0.545) / 2.5;
if (M >= 1.136) and (M <= 5.613) then Cz_bet : = -(-0.056*M + 0.586) / 2.5;
if (M >= 5.613) then Cz_bet : = -0.27924 / 2.5;
end;
if (abs(bet_param) >= 8 * pi / 180) and (abs(bet_param) < 12.5*pi / 180) then begin
	if (M >= 0) and (M <= 0.301) then Cz_bet : = -0.62224 / 2.5;
if (M >= 0.301) and (M <= 0.954) then Cz_bet : = -(0.015*M + 0.618) / 2.5;
if (M >= 0.954) and (M <= 1.136) then Cz_bet : = -(1.157*M - 0.473) / 2.5;
if (M >= 1.136) and (M <= 5.613) then Cz_bet : = -(0.012*M*M*M - 0.121*M*M + 0.302*M - 0.663) / 2.5;
if (M >= 5.613) then Cz_bet : = -0.72451 / 2.5;
end;
if (abs(bet_param) >= 12.5*pi / 180) and (abs(bet_param) <= 17.5*pi / 180) then begin
if (M >= 0) and (M <= 0.301) then Cz_bet : = -0.98213 / 2.5;
if (M >= 0.301) and (M <= 0.954) then Cz_bet : = -(0.164*M + 0.933) / 2.5;
if (M >= 0.954) and (M <= 1.316) then Cz_bet : = -(-2.552*M*M + 7.025*M - 3.291) / 2.5;
if (M >= 1.316) and (M <= 5.613) then Cz_bet : = -(0.002*M*M - 0.154*M + 1.763) / 2.5;
if (M >= 5.613) then Cz_bet : = -0.99124 / 2.5;
end;
if (abs(bet_param) >= 17.5*pi / 180) and (abs(bet_param) <= 25 * pi / 180) then begin
if (M >= 0) and (M <= 0.301) then Cz_bet : = -1.39881 / 2.5;
if (M >= 0.301) and (M <= 0.954) then Cz_bet : = -(1.154*M*M - 1.103*M + 1.626) / 2.5;
if (M >= 0.954) and (M <= 1.316) then Cz_bet : = -(-6.024*M*M + 15.74*M - 7.917) / 2.5;
if (M >= 1.316) and (M <= 5.613) then Cz_bet : = -(0.042*M*M - 0.51*M + 3.042) / 2.5;
if (M >= 5.613) then Cz_bet : = -1.54017 / 2.5;
end;
if (abs(bet_param) >= 25 * pi / 180) and (abs(bet_param) <= 37.5*pi / 180) then begin
if (M >= 0) and (M <= 0.301) then Cz_bet : = -2.33457 / 2.5;
if (M >= 0.301) and (M <= 1.136) then Cz_bet : = -(1.88*M*M - 0.802*M + 2.406) / 2.5;
if (M >= 1.136) and (M <= 1.695) then Cz_bet : = -(-1.778*M*M + 4.936*M + 0.62) / 2.5;
if (M >= 1.695) and (M <= 5.613) then Cz_bet : = -(0.08*M*M - 0.855*M + 5.098) / 2.5;
if (M >= 5.613) then Cz_bet : = -2.83184 / 2.5;
end;
if (abs(bet_param) >= 37.5*pi / 180) and (abs(bet_param) <= 52.5*pi / 180) then begin
if (M >= 0) and (M <= 0.301) then Cz_bet : = -3.79029 / 2.5;
if (M >= 0.301) and (M <= 0.702) then Cz_bet : = -(-0.05*M + 3.805) / 2.5;
if (M >= 0.702) and (M <= 1.316) then Cz_bet : = -(-10.29*M*M + 24.8*M - 8.572) / 2.5;
if (M >= 1.316) and (M <= 5.613) then Cz_bet : = -(-0.02*M*M*M + 0.317*M*M - 1.663*M + 7.939) / 2.5;
if (M >= 5.613) then Cz_bet : = -4.96039 / 2.5;
end;
if (abs(bet_param) >= 52.5*pi / 180) and (abs(bet_param) <= 70 * pi / 180) then begin
if (M >= 0) and (M <= 0.301) then Cz_bet : = -(3.35042) / 2.5;
if (M >= 0.301) and (M <= 0.954) then Cz_bet : = -(-0.166*M*M + 2.101*M + 2.733) / 2.5;
if (M >= 0.954) and (M <= 1.136) then Cz_bet : = -(19.36*M - 13.88) / 2.5;
if (M >= 1.136) and (M <= 1.502) then Cz_bet : = -(-5.248*M*M + 13.12*M - 0.029) / 2.5;
if (M >= 1.502) and (M <= 5.613) then Cz_bet : = -(-0.034*M*M*M + 0.447*M*M - 1.995*M + 9.937) / 2.5;
if (M >= 5.613) then Cz_bet : = -6.80484 / 2.5;
end;
if (abs(bet_param) >= 70 * pi / 180) and (abs(bet_param) <= 90 * pi / 180) then begin
if (M >= 0) and (M <= 0.301) then Cz_bet : = -3.34229 / 2.5;
if (M >= 0.301) and (M <= 0.702) then Cz_bet : = -(2.442*M + 2.607) / 2.5;
if (M >= 0.702) and (M <= 1.136) then Cz_bet : = -(-8.778*M*M + 26.08*M - 9.663) / 2.5;
if (M >= 1.136) and (M <= 3.366) then Cz_bet : = -(-0.123*M*M*M + 0.891*M*M - 2.513*M + 10.56) / 2.5;
if (M >= 3.366) and (M <= 3.995) then Cz_bet : = -(-0.186*M*M + 1.049*M + 6.091) / 2.5;
if (M >= 3.995) and (M <= 5.613) then Cz_bet : = -(0.065*M + 7.144) / 2.5;
if (M >= 5.613) then Cz_bet : = -7.48732 / 2.5;
end;
if Hg > 80000 then Cz_bet : = 0;
end;*/




