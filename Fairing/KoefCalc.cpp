#include "pch.h"
#include "KoefCalc.h"


/////////////////////


double CyAlPas(double mach, double alpha, double Hg) {
	double Cy_alf = 0;

	double alf_Deg = alpha * toDeg;
	double M = mach;

	double Cy1, Cy2;

	if (mach < 0.8) {
		if (alpha*toDeg < -90) { Cy_alf = -1.5733E-5*alf_Deg*alf_Deg*alf_Deg - 5.08024E-3 * alf_Deg *alf_Deg - 0.53082*alf_Deg - 24.31798; 
		}
		else {
			if (alpha*toDeg < 30) {
				Cy_alf = -1.248959E-5*alf_Deg*alf_Deg*alf_Deg - 5.7044E-4 * alf_Deg *alf_Deg + 0.1239*alf_Deg + 0.4279;
			}
			else {
				Cy_alf = -0.33309E-3 * alf_Deg *alf_Deg + 0.036895*alf_Deg + 2.4855;
			}
		}
	} else {
		if (mach < 1.2) {
			if (alpha*toDeg < -90) {
				Cy1 = -1.5733E-5*alf_Deg*alf_Deg*alf_Deg - 5.08024E-3 * alf_Deg *alf_Deg - 0.53082*alf_Deg - 24.31798;
			}
			else {
				if (alpha*toDeg < 30) {
					Cy1 = -1.248959E-5*alf_Deg*alf_Deg*alf_Deg - 5.7044E-4 * alf_Deg *alf_Deg + 0.1239*alf_Deg + 0.4279;
				}
				else {
					Cy1 = -0.33309E-3 * alf_Deg *alf_Deg + 0.036895*alf_Deg + 2.4855;
				}
			}

			if (alpha*toDeg < -80) {
				Cy2 = 1.177E-3 * alf_Deg *alf_Deg + 0.2448*alf_Deg + 4.1653;
			}
			else {
				if (alpha*toDeg < 40) {
					Cy2 = 2.618E-7*alf_Deg*alf_Deg*alf_Deg*alf_Deg +7.493E-6*alf_Deg*alf_Deg*alf_Deg - 8.893E-4 * alf_Deg *alf_Deg + 0.1242*alf_Deg + 0.5675;
				}
				else {
					Cy2 = -0.8E-3 * alf_Deg *alf_Deg + 0.1211*alf_Deg + 1.3993;
				}
			}
			Cy_alf = (Cy2*(mach - 0.8) + Cy1 * (1.2 - mach)) / (1.2 - 0.8);
		}
		else {
			if (mach < 2.0) {
				if (alpha*toDeg < -60) {
					Cy2 = 0.9258E-3 * alf_Deg *alf_Deg + 0.1753*alf_Deg + 0.008778;
				}
				else {
					if (alpha*toDeg < 30) {
						Cy2 = - 4E-4 * alf_Deg *alf_Deg + 0.1126*alf_Deg + 0.6066;
					}
					else {
						Cy2 = -0.7674E-3 * alf_Deg *alf_Deg + 0.1244*alf_Deg + 0.6897;
					}
				}

				if (alpha*toDeg < -80) {
					Cy1 = 1.177E-3 * alf_Deg *alf_Deg + 0.2448*alf_Deg + 4.1653;
				}
				else {
					if (alpha*toDeg < 40) {
						Cy1 = 2.618E-7*alf_Deg*alf_Deg*alf_Deg*alf_Deg + 7.493E-6*alf_Deg*alf_Deg*alf_Deg - 8.893E-4 * alf_Deg *alf_Deg + 0.1242*alf_Deg + 0.5675;
					}
					else {
						Cy1 = -0.8E-3 * alf_Deg *alf_Deg + 0.1211*alf_Deg + 1.3993;
					}
				}
				Cy_alf = (Cy2*(mach - 1.2) + Cy1 * (2.0 - mach)) / (2.0 - 1.2);
			}
			else {
				if (mach < 4.5) {
					if (alpha*toDeg < -60) {
						Cy1 = 0.9258E-3 * alf_Deg *alf_Deg + 0.1753*alf_Deg + 0.008778;
					}
					else {
						if (alpha*toDeg < 30) {
							Cy1 = -4E-4 * alf_Deg *alf_Deg + 0.1126*alf_Deg + 0.6066;
						}
						else {
							Cy1 = -0.7674E-3 * alf_Deg *alf_Deg + 0.1244*alf_Deg + 0.6897;
						}
					}

					if (alpha*toDeg < -90) {
						Cy2 = 1.05E-3 * alf_Deg *alf_Deg + 0.2102*alf_Deg + 2.9173;
					}
					else {
						if (alpha*toDeg < 10) {
							Cy2 = -2.504E-5*alf_Deg*alf_Deg*alf_Deg - 2.748E-3 * alf_Deg *alf_Deg + 0.0486*alf_Deg + 0.8235;
						}
						else {
							Cy2 = 5.6256E-8*alf_Deg*alf_Deg*alf_Deg*alf_Deg - 1.855E-5*alf_Deg*alf_Deg*alf_Deg + 11.46E-4 * alf_Deg *alf_Deg + 0.0613*alf_Deg + 0.3039;
						}
					}
					Cy_alf = (Cy2*(mach - 2.0) + Cy1 * (4.5 - mach)) / (4.5 - 2.0);
				}
				else {
					if (alpha*toDeg < -90) {
						Cy_alf = 1.05E-3 * alf_Deg *alf_Deg + 0.2102*alf_Deg + 2.9173;
					}
					else {
						if (alpha*toDeg < 10) {
							Cy_alf = -2.504E-5*alf_Deg*alf_Deg*alf_Deg - 2.748E-3 * alf_Deg *alf_Deg + 0.0486*alf_Deg + 0.8235;
						}
						else {
							Cy_alf = 5.6256E-8*alf_Deg*alf_Deg*alf_Deg*alf_Deg - 1.855E-5*alf_Deg*alf_Deg*alf_Deg + 11.46E-4 * alf_Deg *alf_Deg + 0.0613*alf_Deg + 0.3039;
						}
					}
				}
			}

		}


		
	}
	

	
	return Cy_alf;
}





double CzBettaPas(double mach, double alpha, double betta) {
	double Cz = 0;
	double Cz1 = 0;
	double Cz2 = 0;

	double alphaDeg = alpha * toDeg;
	double alphaDeg2 = alphaDeg * alphaDeg;
	double bettaA = abs(betta);
	double betta1 = (bettaA > PI / 2) ? PI - bettaA : bettaA;
	double bettaAbsDeg = abs(betta1 * toDeg);
	if (bettaAbsDeg < 30) {
		Cz1 = 0;
		if (alphaDeg < 130.485) {
			Cz2 = 0.20918679 - 0.044540495*alphaDeg + 2.9081254E-4*alphaDeg2;
		}
		else {
			if (alpha * toDeg < 200) {
				Cz2 = -7.071903494 + 0.07632497*alphaDeg - 2.04505654E-4*alphaDeg2;
			}
			else {
				Cz2 = -4.4735130585 + 0.034440515*alphaDeg - 6.1251803E-5*alphaDeg2;
			}
		}

		Cz = (Cz2*(bettaAbsDeg - 0.0) + Cz1 * (30 - bettaAbsDeg)) / (30 - 0.0);
		Cz = (betta < 0) ? -Cz : Cz;
	}
	else {
		if (abs(betta1)*toDeg < 60) {
			if (alphaDeg < 130.485) {
				Cz1 = 0.20918679 - 0.044540495*alphaDeg + 2.9081254E-4*alphaDeg2;
			}
			else {
				if (alphaDeg < 200) {
					Cz1 = -7.071903494 + 0.07632497*alphaDeg - 2.04505654E-4*alphaDeg2;
				}
				else {
					Cz1 = -4.4735130585 + 0.034440515*alphaDeg - 6.1251803E-5*alphaDeg2;
				}
			}
			if (alphaDeg < 131.068) {
				Cz2 = 0.470694597 - 0.07891926*alphaDeg + 4.945039E-4*alphaDeg2;
			}
			else {
				if (alphaDeg < 200) {
					Cz2 = -13.35566025 + 0.1419816*alphaDeg - 3.7534874E-4*alphaDeg2;
				}
				else {
					Cz2 = -19.49398763 + 0.151225153*alphaDeg - 2.713184E-4*alphaDeg2;
				}
			}
			Cz = (Cz2*(bettaAbsDeg - 30.0) + Cz1 * (60 - bettaAbsDeg)) / (60 - 30.0);
			Cz = (betta < 0) ? -Cz : Cz;
		}else{
			if (abs(betta1)*toDeg < 90) {
				if (alphaDeg < 131.068) {
					Cz1 = 0.470694597 - 0.07891926*alphaDeg + 4.945039E-4*alphaDeg2;
				}
				else {
					if (alphaDeg < 200) {
						Cz1 = -13.35566025 + 0.1419816*alphaDeg - 3.7534874E-4*alphaDeg2;
					}
					else {
						Cz1 = -19.49398763 + 0.151225153*alphaDeg - 2.713184E-4*alphaDeg2;
					}
				}


				if (alphaDeg < 130.485) {
					Cz2 = 0.670436646 - 0.08099633*alphaDeg + 5.05206E-4*alphaDeg2;
				}
				else {
					if (alphaDeg < 200) {
						Cz2 = -12.556949 + 0.13115159*alphaDeg - 3.384524E-4*alphaDeg2;
					}
					else {
						Cz2 = -31.49609359 + 0.24447658*alphaDeg - 4.3896712E-4*alphaDeg2;
					}
				}
				Cz = (Cz2*(bettaAbsDeg - 60.0) + Cz1 * (90 - bettaAbsDeg)) / (90 - 60.0);
				Cz = (betta < 0) ? -Cz : Cz;

			}
		}

		
	}

	return Cz;
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

	double alf_Deg = alpha * toDeg;
	double M = mach;

	double Cx1, Cx2;

	if (mach < 0.8) {
		if (alpha*toDeg < -90) {
			fCx = -5.399905E-6*alf_Deg*alf_Deg*alf_Deg - 0.002028125 * alf_Deg *alf_Deg - 0.2471424*alf_Deg - 11.002250685998;
		}
		else {
			if (alpha*toDeg < 0) {
				fCx = 3.7729423E-4 * alf_Deg *alf_Deg + 0.055806*alf_Deg + 0.69685568;
			}
			else {
				fCx = 4.01243E-7 * alf_Deg *alf_Deg*alf_Deg - 1.4709E-4 * alf_Deg *alf_Deg + 0.0056699*alf_Deg + 0.6833;
			}
		}
	}
	else {
		if (mach < 1.2) {
			if (alpha*toDeg < -90) {
				Cx1 = -5.399905E-6*alf_Deg*alf_Deg*alf_Deg - 0.002028125 * alf_Deg *alf_Deg - 0.2471424*alf_Deg - 11.002250685998;
			}
			else {
				if (alpha*toDeg < 0) {
					Cx1 = 3.7729423E-4 * alf_Deg *alf_Deg + 0.055806*alf_Deg + 0.69685568;
				}
				else {
					Cx1 = 4.01243E-7 * alf_Deg *alf_Deg*alf_Deg - 1.4709E-4 * alf_Deg *alf_Deg + 0.0056699*alf_Deg + 0.6833;
				}
			}

			if (alpha*toDeg < -120) {
				Cx2 = 0.4563E-3 * alf_Deg *alf_Deg + 0.128*alf_Deg + 7.3263;
			}
			else {
				if (alpha*toDeg < 0) {
					Cx2 = -1.14E-7*alf_Deg*alf_Deg*alf_Deg*alf_Deg - 2.8189E-5*alf_Deg*alf_Deg*alf_Deg - 1.8595E-3 * alf_Deg *alf_Deg + 0.00463*alf_Deg + 0.7628;
				}
				else {
					Cx2 = -0.1267E-3 * alf_Deg *alf_Deg + 0.0114*alf_Deg + 0.9322;
				}
			}
			fCx = (Cx2*(mach - 0.8) + Cx1 * (1.2 - mach)) / (1.2 - 0.8);
		}
		else {
			if (mach < 2.0) {
				if (alpha*toDeg < -110) {
					Cx2 = 0.23458E-3 * alf_Deg *alf_Deg + 0.0605*alf_Deg + 2.3739;
				}
				else {
					if (alpha*toDeg < 10) {
						Cx2 = -3.184E-8*alf_Deg*alf_Deg*alf_Deg*alf_Deg - 1.03E-5*alf_Deg*alf_Deg*alf_Deg - 0.7903E-3 * alf_Deg *alf_Deg + 0.01432*alf_Deg + 0.6342;
					}
					else {
						Cx2 = -0.126E-3 * alf_Deg *alf_Deg + 0.0129*alf_Deg + 0.7442;
					}
				}

				if (alpha*toDeg < -120) {
					Cx1 = 0.4563E-3 * alf_Deg *alf_Deg + 0.128*alf_Deg + 7.3263;
				}
				else {
					if (alpha*toDeg < 0) {
						Cx1 = -1.14E-7*alf_Deg*alf_Deg*alf_Deg*alf_Deg - 2.8189E-5*alf_Deg*alf_Deg*alf_Deg - 1.8595E-3 * alf_Deg *alf_Deg + 0.00463*alf_Deg + 0.7628;
					}
					else {
						Cx1 = -0.1267E-3 * alf_Deg *alf_Deg + 0.0114*alf_Deg + 0.9322;
					}
				}
				fCx = (Cx2*(mach - 1.2) + Cx1 * (2.0 - mach)) / (2.0 - 1.2);
			}
			else {
				if (mach < 4.5) {
					if (alpha*toDeg < -110) {
						Cx1 = 0.23458E-3 * alf_Deg *alf_Deg + 0.0605*alf_Deg + 2.3739;
					}
					else {
						if (alpha*toDeg < 10) {
							Cx1 = -3.184E-8*alf_Deg*alf_Deg*alf_Deg*alf_Deg - 1.03E-5*alf_Deg*alf_Deg*alf_Deg - 0.7903E-3 * alf_Deg *alf_Deg + 0.01432*alf_Deg + 0.6342;
						}
						else {
							Cx1 = -0.126E-3 * alf_Deg *alf_Deg + 0.0129*alf_Deg + 0.7442;
						}
					}

					if (alpha*toDeg < -70) {
						Cx2 = 1.597586E-7*alf_Deg*alf_Deg*alf_Deg*alf_Deg + 7.48421E-5*alf_Deg*alf_Deg*alf_Deg + 0.01275 * alf_Deg *alf_Deg + 0.936455*alf_Deg + 23.5509;
					}
					else {
						if (alpha*toDeg < 0) {
							Cx2 = 2.1356E-6*alf_Deg*alf_Deg*alf_Deg - 1.996E-4 * alf_Deg *alf_Deg + 0.00024*alf_Deg + 0.2909;
						}
						else {
							Cx2 = 9.381E-7*alf_Deg*alf_Deg*alf_Deg - 3.7364E-4 * alf_Deg *alf_Deg + 0.0321*alf_Deg + 0.2136;
						}
					}
					fCx = (Cx2*(mach - 2.0) + Cx1 * (4.5 - mach)) / (4.5 - 2.0);
				}
				else {
					if (alpha*toDeg < -70) {
						fCx = 1.597586E-7*alf_Deg*alf_Deg*alf_Deg*alf_Deg + 7.48421E-5*alf_Deg*alf_Deg*alf_Deg + 0.01275 * alf_Deg *alf_Deg + 0.936455*alf_Deg + 23.5509;
					}
					else {
						if (alpha*toDeg < 0) {
							fCx = 2.1356E-6*alf_Deg*alf_Deg*alf_Deg - 1.996E-4 * alf_Deg *alf_Deg + 0.00024*alf_Deg + 0.2909;
						}
						else {
							fCx = 9.381E-7*alf_Deg*alf_Deg*alf_Deg - 3.7364E-4 * alf_Deg *alf_Deg + 0.0321*alf_Deg + 0.2136;
						}
					}
				}
			}

		}



	}
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
	double alphaDeg = alpha * toDeg;
	
	if (alphaDeg < 160) {
		Mz_wz = -1.297788E-8 * alphaDeg * alphaDeg * alphaDeg + 3.722716E-5 * alphaDeg * alphaDeg - 0.00580598 *alphaDeg - 0.0191;
	}
	else {
		if (alphaDeg < 180) {
			Mz_wz = -0.0014244 *alphaDeg + 0.16869;
		}
		else {
			double a2 = alphaDeg*alphaDeg;
			double a4 = a2 * a2;
			Mz_wz = -11.048193 + 0.20715736*alphaDeg - 0.00136936488*a2 + 3.74171187E-6*a2*alphaDeg - 3.61118363E-9*a4;
		}
	}

	
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

	double alphaDeg = alpha * toDeg;
	double M = mach;

	double mz1, mz2;

	
	double alphaDeg2 = alphaDeg * alphaDeg;
	double alphaDeg3 = alphaDeg2 * alphaDeg;
	double alphaDeg4 = alphaDeg2 * alphaDeg2;
	double alphaDeg5 = alphaDeg3 * alphaDeg2;

	if (mach < 0.8) {
		if (alpha*toDeg < -90) {
			Mz_alf = -1.484E-4 * alphaDeg2 - 0.0413*alphaDeg - 2.6698;
		}
		else {
			if (alpha*toDeg < 60) {
				Mz_alf = 1.52666E-8*alphaDeg4 - 1.04736E-6*alphaDeg3 - 1.67E-4 * alphaDeg2 + 0.0093*alphaDeg + 0.2788;
			}
			else {
				Mz_alf = 5.924E-7 * alphaDeg3 - 1.021E-4 * alphaDeg2 - 0.005386*alphaDeg + 0.786;
			}
		}
	}
	else {
		if (mach < 1.2) {
			if (alpha*toDeg < -90) {
				mz1 = -1.484E-4 * alphaDeg2 - 0.0413*alphaDeg - 2.6698;
			}
			else {
				if (alpha*toDeg < 60) {
					mz1 = 1.52666E-8*alphaDeg4 - 1.04736E-6*alphaDeg3 - 1.67E-4 * alphaDeg2 + 0.0093*alphaDeg + 0.2788;
				}
				else {
					mz1 = 5.924E-7 * alphaDeg3 - 1.021E-4 * alphaDeg2 - 0.005386*alphaDeg + 0.786;
				}
			}

			if (alpha*toDeg < -70) {
				mz2 = -0.13303E-3 * alphaDeg2 - 0.0372*alphaDeg - 2.42997;
			}
			else {
				if (alpha*toDeg < 50) {
					mz2 = -1.4721E-6*alphaDeg3 - 0.13064E-3 * alphaDeg2 + 0.009851*alphaDeg + 0.36698;
				}
				else {
					mz2 = 1.873E-6*alphaDeg3 - 0.5487E-3 * alphaDeg2 + 0.0408212*alphaDeg - 0.56743;
				}
			}
			Mz_alf = (mz2*(mach - 0.8) + mz1 * (1.2 - mach)) / (1.2 - 0.8);
		}
		else {
			if (mach < 2.0) {
				if (alpha*toDeg < -80) {
					mz2 = -0.14842E-3 * alphaDeg2 - 0.040067*alphaDeg - 2.5726;
				}
				else {
					if (alpha*toDeg < 70) {
						mz2 = -1.5795E-10*alphaDeg5 + 3.2E-9*alphaDeg4 + 6.846E-7*alphaDeg3 - 8.737E-5 * alphaDeg2 + 0.00407*alphaDeg + 0.2857;
					}
					else {
						mz2 = 0.107E-3 * alphaDeg2 - 0.03029*alphaDeg + 1.8446;
					}
				}

				if (alpha*toDeg < -70) {
					mz1 = -0.13303E-3 * alphaDeg2 - 0.0372*alphaDeg - 2.42997;
				}
				else {
					if (alpha*toDeg < 50) {
						mz1 = -1.4721E-6*alphaDeg3 - 0.13064E-3 * alphaDeg2 + 0.009851*alphaDeg + 0.36698;
					}
					else {
						mz1 = 1.873E-6*alphaDeg3 - 0.5487E-3 * alphaDeg2 + 0.0408212*alphaDeg - 0.56743;
					}
				}
				Mz_alf = (mz2*(mach - 1.2) + mz1 * (2.0 - mach)) / (2.0 - 1.2);
			}
			else {
				if (mach < 4.5) {
					if (alpha*toDeg < -80) {
						mz1 = -0.14842E-3 * alphaDeg2 - 0.040067*alphaDeg - 2.5726;
					}
					else {
						if (alpha*toDeg < 70) {
							mz1 = -1.5795E-10*alphaDeg5 +3.2E-9*alphaDeg4 + 6.846E-7*alphaDeg3 - 8.737E-5 * alphaDeg2 + 0.00407*alphaDeg + 0.2857;
						}
						else {
							mz1 = 0.107E-3 * alphaDeg2 - 0.03029*alphaDeg + 1.8446;
						}
					}

					if (alpha*toDeg < -100) {
						mz2 = 4.04756E-7*alphaDeg4 + 2.1536E-4*alphaDeg3 + 0.04208 * alphaDeg2 + 3.58149*alphaDeg + 112.245;
					}
					else {
						if (alpha*toDeg < 0) {
						//	mz2 = 6.6703E-8*alf_Deg*alf_Deg*alf_Deg*alf_Deg + 9.1283E-6*alf_Deg*alf_Deg*alf_Deg + 2.428E-4 * alf_Deg *alf_Deg + 6.3847E-4*alf_Deg + 0.12766;
							mz2 = -1.6063E-8*alphaDeg4 - 6.325E-6*alphaDeg3 - 6.432E-4 * alphaDeg2 + 1.514E-2*alphaDeg + 0.0921;
						}
						else {
							mz2 = -5.0335E-9*alphaDeg4 + 2.5308E-6*alphaDeg3 - 0.000398 * alphaDeg2 + 0.0178*alphaDeg + 0.1006;
						}
					}
					Mz_alf = (mz2*(mach - 2.0) + mz1 * (4.5 - mach)) / (4.5 - 2.0);
				} else {
					if (alpha*toDeg < -100) {
						Mz_alf = 4.04756E-7*alphaDeg4 + 2.1536E-4*alphaDeg3 + 0.04208 * alphaDeg2 + 3.58149*alphaDeg + 112.245;
					}
					else {
						if (alpha*toDeg < 0) {
							Mz_alf = -1.6063E-8*alphaDeg4 - 6.325E-6*alphaDeg3 - 6.432E-4 * alphaDeg2 - 1.514E-2*alphaDeg + 0.0921;
						}
						else {
							Mz_alf = -5.0335E-9*alphaDeg4 + 2.5308E-6*alphaDeg3 - 0.000398 * alphaDeg2 + 0.0178*alphaDeg + 0.1006;
						}
					}
				}
			}

		}



	}
	return Mz_alf;
}




double MyOmegaY(double mach, double betta) {
	return 0;

}



double MyBettaPas(double mach, double alpha, double betta) {

	double my = 0;
	double my1 = 0;
	double my2 = 0;
	double bettaA = abs(betta);
	double betta1 = (bettaA > PI / 2) ? PI - bettaA : bettaA;
	double alphaDeg = alpha * toDeg;
	double bettaAbsDeg = abs(betta1*toDeg);
	double alphaDeg2 = alphaDeg * alphaDeg;
	double alphaDeg3 = alphaDeg2 * alphaDeg;
	double alphaDeg4 = alphaDeg2 * alphaDeg2;
	double alphaDeg5 = alphaDeg3 * alphaDeg2;

	if (bettaAbsDeg < 30) {
		my1 = 0;
		if (alpha * toDeg < 181.25) {
			my2 = 1.177265E-3 + 2.0660087E-3*alphaDeg + 5.1360683E-5*alphaDeg2 - 1.99451668E-6*alphaDeg3 + 1.61694008413613E-8*alphaDeg4 - 3.90498782894E-11*alphaDeg5;
		}
		else {
			
				my2 = -6.59739187 + 0.13780116*alphaDeg - 1.151492833125E-3*alphaDeg2 + 4.801379052E-6*alphaDeg3 - 9.9470177003001E-9*alphaDeg4 + 8.15104025386042E-12*alphaDeg5;
			
		}

		my = (my2*(bettaAbsDeg - 0.0) + my1 * (30 - bettaAbsDeg)) / (30 - 0.0);
		my = (betta < 0) ? -my : my;
	}
	else {
		if (bettaAbsDeg < 60) {
			if (alpha * toDeg < 181.25) {
				my1 = 1.177265E-3 + 2.0660087E-3*alphaDeg + 5.1360683E-5*alphaDeg2 - 1.99451668E-6*alphaDeg3 + 1.61694008413613E-8*alphaDeg4 - 3.90498782894E-11*alphaDeg5;
			}
			else {				
					my1 = -6.59739187 + 0.13780116*alphaDeg - 1.151492833125E-3*alphaDeg2 + 4.801379052E-6*alphaDeg3 - 9.9470177003001E-9*alphaDeg4 + 8.15104025386042E-12*alphaDeg5;				
			}
			if (alpha * toDeg < 161.818) {
				my2 = 1.293096395E-3 + 3.4148346676E-3*alphaDeg + 4.7291469E-5*alphaDeg2- 2.4200543938E-6*alphaDeg3 + 2.017695332E-8*alphaDeg4 - 4.85801728813182E-11*alphaDeg5;
			}
			else {				
					my2 = 6.0075128 - 0.1103382468*alphaDeg + 7.33545E-4*alphaDeg2 - 2.089793715E-6*alphaDeg3 + 2.171407611E-9*alphaDeg4 - 5.2208322772E-14*alphaDeg5;				
			}
			my = (my2*(bettaAbsDeg - 30.0) + my1 * (60 - bettaAbsDeg)) / (60 - 30.0);
			my = (betta < 0) ? -my : my;
		}
		else {
			if (bettaAbsDeg < 90) {
				if (alpha * toDeg < 161.818) {
					my1 = 1.293096395E-3 + 3.4148346676E-3*alphaDeg + 4.7291469E-5*alphaDeg2 - 2.4200543938E-6*alphaDeg3 + 2.017695332E-8*alphaDeg4 - 4.85801728813182E-11*alphaDeg5;
				}
				else {
					
						my1 = 6.0075128 - 0.1103382468*alphaDeg + 7.33545E-4*alphaDeg2 - 2.089793715E-6*alphaDeg3 + 2.171407611E-9*alphaDeg4 - 5.2208322772E-14*alphaDeg5;					
				}


				if (alpha * toDeg < 180.625) {
					my2 = -4.1483153072E-3 + 5.32512481E-3*alphaDeg - 4.25749351E-5*alphaDeg2 - 1.0511163881E-6*alphaDeg3 + 1.2086021876E-8*alphaDeg4 - 3.24614262042153E-11*alphaDeg5;
				}
				else {
					
						my2 = 50.089004333 - 0.9840487*alphaDeg + 7.519415035E-3*alphaDeg2 - 2.7900414507E-5*alphaDeg3 + 5.02879868E-8*alphaDeg4 - 3.527056354051E-11*alphaDeg5;					
				}
				my = (my2*(bettaAbsDeg - 60.0) + my1 * (90 - bettaAbsDeg)) / (90 - 60.0);
				my = (betta < 0) ? -my : my;

			}
		}
	}


	return my;
}

double MxBettaPas(double mach, double alpha, double betta) {

	double mx = 0;
	double mx1 = 0;
	double mx2 = 0;
	double bettaA = abs(betta);
	double betta1 = (bettaA > PI / 2) ? PI - bettaA : bettaA;
	double alphaDeg = alpha * toDeg;
	double bettaAbsDeg = abs(betta1*toDeg);
	double alphaDeg2 = alphaDeg * alphaDeg;
	double alphaDeg3 = alphaDeg2 * alphaDeg;
	double alphaDeg4 = alphaDeg2 * alphaDeg2;
	double alphaDeg5 = alphaDeg3 * alphaDeg2;
	if (bettaAbsDeg < 30) {
		mx1 = 0;
		if (alpha * toDeg < 180) {
			mx2 = -4.11595338E-3 + 4.322027523E-4*alphaDeg - 1.518620114E-4*alphaDeg2 + 2.3259526793E-6*alphaDeg3 - 1.2024717703E-8*alphaDeg4 + 2.06506966141321E-11*alphaDeg5;
		}
		else {

			mx2 = 6.873320132 - 0.1249343665*alphaDeg + 8.830828935E-4*alphaDeg2 - 3.04725145156E-6*alphaDeg3 + 5.17817218653E-9*alphaDeg4 - 3.49735022619433E-12*alphaDeg5;

		}

		mx = (mx2*(bettaAbsDeg - 0.0) + mx1 * (30 - bettaAbsDeg)) / (30 - 0.0);
		mx = (betta < 0) ? -mx : mx;
	}
	else {
		if (bettaAbsDeg < 60) {
			if (alpha * toDeg < 180) {
				mx1 = -4.11595338E-3 + 4.322027523E-4*alphaDeg - 1.518620114E-4*alphaDeg2 + 2.3259526793E-6*alphaDeg3 - 1.2024717703E-8*alphaDeg4 + 2.06506966141321E-11*alphaDeg5;
			}
			else {
				mx1 = 6.873320132 - 0.1249343665*alphaDeg + 8.830828935E-4*alphaDeg2 - 3.04725145156E-6*alphaDeg3 + 5.17817218653E-9*alphaDeg4 - 3.49735022619433E-12*alphaDeg5;
			}
			if (alpha * toDeg < 178.812) {
				mx2 = -4.819184741E-3 + 3.720456549E-4*alphaDeg - 2.3052581254E-4*alphaDeg2 + 3.29557444863E-6*alphaDeg3 - 1.529297242E-8*alphaDeg4 + 2.24075890939032E-11*alphaDeg5;
			}
			else {

				mx2 = 13.56323183 - 0.22266758673*alphaDeg + 1.34401620049E-3*alphaDeg2 - 3.61050893514E-6*alphaDeg3 + 4.0458625237E-9*alphaDeg4 - 1.1729417603405E-12*alphaDeg5;

			}
			mx = (mx2*(bettaAbsDeg - 30.0) + mx1 * (60 - bettaAbsDeg)) / (60 - 30.0);
			mx = (betta < 0) ? -mx : mx;
		}
		else {
			if (bettaAbsDeg < 90) {
				if (alpha * toDeg < 178.812) {
					mx1 = -4.819184741E-3 + 3.720456549E-4*alphaDeg - 2.3052581254E-4*alphaDeg2 + 3.29557444863E-6*alphaDeg3 - 1.529297242E-8*alphaDeg4 + 2.24075890939032E-11*alphaDeg5;
				}
				else {

					mx1 = 13.56323183 - 0.22266758673*alphaDeg + 1.34401620049E-3*alphaDeg2 - 3.61050893514E-6*alphaDeg3 + 4.0458625237E-9*alphaDeg4 - 1.1729417603405E-12*alphaDeg5;

				}


				if (alpha * toDeg < 180) {
					mx2 = -4.9800813517E-3 - 3.76518127E-4*alphaDeg - 1.862149644E-4*alphaDeg2 + 2.6486065044E-6*alphaDeg3 - 1.159565965E-8*alphaDeg4 + 1.493840398812913E-11*alphaDeg5;
				}
				else {

					mx2 = 9.135795554 - 0.10251782329*alphaDeg + 1.330991724E-4*alphaDeg2 + 2.1199874904E-6*alphaDeg3 - 8.7538390454E-9*alphaDeg4 + 9.69848286137204E-12*alphaDeg5;

				}
				mx = (mx2*(bettaAbsDeg - 60.0) + mx1 * (90 - bettaAbsDeg)) / (90 - 60.0);
				mx = (betta < 0) ? -mx : mx;

			}
		}
	}


	return mx;
}
