#include "pch.h"
#include "KoefCalc.h"


/////////////////////


double CyAlPas(double mach, double alpha, double Hg) {
	double Cy = 0;

	double alf_Deg = alpha * toDeg;
	double M = mach;
	double a1 = alpha;
	double a2 = alpha * alpha;
	double a3 = a2 * a1;
	double a4 = a3 * a1;
	double a5 = a4 * a1;
	double a6 = a5 * a1;
	double Cy1, Cy2;

	if (mach < 0.8) {
		if (alpha < 0) { 
			Cy2 = 0.190594304652482 - 6.38406136886791*a1 - 33.5164265635582*a2 - 30.704636866701*a3 - 10.9159604526799*a4 - 1.3746778158162*a5;
			Cy1 = 0.2*Cy2;
		}
		else {
			Cy2 = 0.413481820060867 + 4.71103549294835*a1 - 1.76350613168867*a2;
			Cy1 = 0.2*Cy2;
		}
		Cy = (Cy2*(mach - 0.0) + Cy1 * (0.8 - mach)) / (0.8 - 0.0);
		//Cy = Cy2;
	} else {
		if (mach < 1.2) {
			if (alpha < 0) {
				Cy1 = 0.190594304652482 - 6.38406136886791*a1 - 33.5164265635582*a2 - 30.704636866701*a3 - 10.9159604526799*a4 - 1.3746778158162*a5;
				Cy2 = 0.453857757764376 - 5.43687824483694*a1 - 33.8838359850678*a2 - 30.2119328212675*a3 - 10.4390863455876*a4 - 1.29469660098323*a5;
			}
			else {
				Cy1 = 0.413481820060867 + 4.71103549294835*a1 - 1.76350613168867*a2;
				Cy2 = 0.120082691654333 + 2.50582043064633*a1 + 21.3350685506266*a2 - 35.110871622988*a3 + 22.7455980952535*a4 - 6.75651693112704*a5 + 0.748903857717365*a6;
			}

			Cy = (Cy2*(mach - 0.8) + Cy1 * (1.2 - mach)) / (1.2 - 0.8);
		}
		else {
			if (mach < 2.0) {
				if (alpha < 0) {
					Cy1 = 0.453857757764376 - 5.43687824483694*a1 - 33.8838359850678*a2 - 30.2119328212675*a3 - 10.4390863455876*a4 - 1.29469660098323*a5;
					Cy2 = 0.740036659580577 + 5.76397325556497*a1 - 8.23184438631578*a2 - 9.35855376619993*a3 - 3.20462594417264*a4 - 0.389311584226631*a5;
				}
				else {
					Cy1 = 0.120082691654333 + 2.50582043064633*a1 + 21.3350685506266*a2 - 35.110871622988*a3 + 22.7455980952535*a4 - 6.75651693112704*a5 + 0.748903857717365*a6;
					Cy2 = 0.441581999540912 + 2.46220116703884*a1 + 9.29525221217141*a2 - 10.70759758376*a3 + 4.55991812468966*a4 - 0.914541243489724*a5 + 0.0679474605594143*a6;
				}

				Cy = (Cy2*(mach - 1.2) + Cy1 * (2.0 - mach)) / (2.0 - 1.2);
				//Cy = Cy2;
			}
			else {
				if (mach < 4.5) {
					if (alpha < 0) {
						Cy1 = 0.740036659580577 + 5.76397325556497*a1 - 8.23184438631578*a2 - 9.35855376619993*a3 - 3.20462594417264*a4 - 0.389311584226631*a5;
						Cy2 = 0.603100827419247 - 2.10363901415552*a1 - 20.3536891454537*a2 - 17.4794379073754*a3 - 5.6926213325844*a4 - 0.672518772908521*a5;
					}
					else {
						Cy1 = 0.441581999540912 + 2.46220116703884*a1 + 9.29525221217141*a2 - 10.70759758376*a3 + 4.55991812468966*a4 - 0.914541243489724*a5 + 0.0679474605594143*a6;
						Cy2 = 0.581067565888243 + 0.0365097869043904*a1 + 10.1377522720521*a2 - 8.78588520201121*a3 + 2.58510324922671*a4 - 0.265554377627029*a5;
					}
					Cy = (Cy2*(mach - 2.0) + Cy1 * (4.5 - mach)) / (4.5 - 2.0);
				}
				else {
					if (alpha < 0) {
						Cy = 0.603100827419247 - 2.10363901415552*a1 - 20.3536891454537*a2 - 17.4794379073754*a3 - 5.6926213325844*a4 - 0.672518772908521*a5;
					}
					else {
						Cy = 0.581067565888243 + 0.0365097869043904*a1 + 10.1377522720521*a2 - 8.78588520201121*a3 + 2.58510324922671*a4 - 0.265554377627029*a5;
					}
				}
			}
		}		
	}
	
	
		

	return Cy;
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

	if (mach < 0.8) {
		double Cz1 = 0.2*0.94*Cz;
		double Cz2 = 0.94*Cz;
		Cz = (Cz2*(mach - 0.0) + Cz1 * (0.8 - mach)) / (0.8 - 0.0);
		//Cz = Cz2;
	}
	else {
		if (mach < 1.2) {
			double Cz1 = 0.94*Cz;
			double Cz2 = 1.47*Cz;
			Cz = (Cz2*(mach - 0.8) + Cz1 * (1.2 - mach)) / (1.2 - 0.8);
		}
		else {
			if (mach < 2.0) {
				double Cz1 = 1.47*Cz;
				double Cz2 = 1.36*Cz;
				Cz = (Cz2*(mach - 1.2) + Cz1 * (2.0 - mach)) / (2.0 - 1.2);

			
			}
			else {
				if (mach < 4.5) {
					double Cz1 = 1.36*Cz;
					double Cz2 = Cz;
					Cz = (Cz2*(mach - 2.0) + Cz1 * (4.5 - mach)) / (4.5 - 2.0);
					
				}
				
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
	double Cx = 0;

	double alf_Deg = alpha * toDeg;
	double M = mach;
	double a1 = alpha;
	double a2 = alpha * alpha;
	double a3 = a2 * a1;
	double a4 = a3 * a1;
	double a5 = a4 * a1;
	double a6 = a5 * a1;
	double Cx1, Cx2;
	////////////////////////////////////////////////////////////////////////////////////////не исправляла......

	if (mach < 0.8) {
		if (alpha < 0) {
			Cx2 = 0.63926001519838 - 0.162277253027975*a1 - 7.77570070224025*a2 - 9.65679071603829*a3 - 4.97317651231006*a4 - 1.18095024713667*a5 - 0.105813731416767*a6;
			Cx1 = 0.5*Cx2;
		}
		else {
			Cx2 = 0.672129227883525 + 0.108693252293213*a1 + 0.906856987527982*a2 - 1.45582628050291*a3 + 0.633634965318283*a4 - 0.0900901742641017*a5;
			Cx1 = Cx2 * 0.5;
		}
		Cx = (Cx2*(mach - 0.0) + Cx1 * (0.8 - mach)) / (0.8 - 0.0);
		//Cx = Cx2;
	}
	else {
		if (mach < 1.2) {
			if (alpha < 0) {
				Cx1 = 0.63926001519838 - 0.162277253027975*a1 - 7.77570070224025*a2 - 9.65679071603829*a3 - 4.97317651231006*a4 - 1.18095024713667*a5 - 0.105813731416767*a6;
				Cx2 = 0.885119321791345 + 2.78160620119536*a1 + 1.67630437725316*a2 + 2.12735144922785*a3 + 1.74782155319269*a4 + 0.61601427287107*a5 + 0.077448880142836*a6;
			}
			else {
				Cx1 = 0.672129227883525 + 0.108693252293213*a1 + 0.906856987527982*a2 - 1.45582628050291*a3 + 0.633634965318283*a4 - 0.0900901742641017*a5;
				Cx2 = 0.921902946026878 + 1.17222088376394*a1 - 1.17219029325695*a2 + 0.34559204364362*a3 - 0.0382404866079392*a4 - 0.0034203187589978*a5;
			}

			Cx = (Cx2*(mach - 0.8) + Cx1 * (1.2 - mach)) / (1.2 - 0.8);
		}
		else {
			if (mach < 2.0) {
				if (alpha < 0) {
					Cx1 = 0.885119321791345 + 2.78160620119536*a1 + 1.67630437725316*a2 + 2.12735144922785*a3 + 1.74782155319269*a4 + 0.61601427287107*a5 + 0.077448880142836*a6;
					Cx2 = 0.365580708328933 + 0.597715976382403*a1 + 0.557140991020652*a2 + 2.52599088374278*a3 + 2.2059056136048*a4 + 0.729553688770536*a5 + 0.0855518761652196*a6;
				}
				else {
					Cx2 = 0.38349254408567 + 0.529713035341652*a1 + 1.2050290101389*a2 - 1.45954297097336*a3 + 0.466491292961852*a4 - 0.0486311373721188*a5;
					Cx1 = 0.921902946026878 + 1.17222088376394*a1 - 1.17219029325695*a2 + 0.34559204364362*a3 - 0.0382404866079392*a4 - 0.0034203187589978*a5;
				}

				Cx = (Cx2*(mach - 1.2) + Cx1 * (2.0 - mach)) / (2.0 - 1.2);
				//Cx = Cx2;
			}
			else {
				if (mach < 4.5) {
					if (alpha < 0) {
						Cx1 = 0.365580708328933 + 0.597715976382403*a1 + 0.557140991020652*a2 + 2.52599088374278*a3 + 2.2059056136048*a4 + 0.729553688770536*a5 + 0.0855518761652196*a6;
						Cx2 = 0.23674654854381 + 0.443779473147453*a1 + 1.0212358993148*a2 + 3.63436948689837*a3 + 2.97046462631658*a4 + 0.934575727572502*a5 + 0.10378560859549*a6;
					}
					else {
						Cx1 = 0.38349254408567 + 0.529713035341652*a1 + 1.2050290101389*a2 - 1.45954297097336*a3 + 0.466491292961852*a4 - 0.0486311373721188*a5;
						Cx2 = 0.253975793231017 + 0.347950469560399*a1 + 2.10148499591583*a2 - 2.50623631387271*a3 + 0.910987463592628*a4 - 0.110539834971781*a5;
					}
					Cx = (Cx2*(mach - 2.0) + Cx1 * (4.5 - mach)) / (4.5 - 2.0);
				}
				else {
					if (alpha < 0) {
						Cx = 0.23674654854381 + 0.443779473147453*a1 + 1.0212358993148*a2 + 3.63436948689837*a3 + 2.97046462631658*a4 + 0.934575727572502*a5 + 0.10378560859549*a6;
					}
					else {
						Cx = 0.253975793231017 + 0.347950469560399*a1 + 2.10148499591583*a2 - 2.50623631387271*a3 + 0.910987463592628*a4 - 0.110539834971781*a5;
					}
				}
			}

		}
	}


	return Cx;

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
	
	double mz = 0;
	

	double mz1, mz2;

	
	double a1 = alpha;
	double a2 = a1 * a1;
	double a3 = a1 * a2;
	double a4 = a1 * a3;
	double a5 = a1 * a4;
	double a6 = a1 * a5;
	double a7 = a1 * a6;
	if (mach < 0.8) {
		if (alpha < 0) {
			mz2 = 0.30779142651144 - 0.277854038854256*a1 - 6.65707437924101*a2 - 10.3248375924987*a3 - 6.78103257929276*a4 - 2.30964522348498*a5 - 0.401738062965339*a6 - 0.0279362656081643*a7;
			mz1 = 0.2*mz2;
		}
		else {
			mz2 = 0.344548567530313 + 0.562705635856*a1 - 0.851458363904556*a2 + 0.610883020728777*a3 - 0.5142261300719*a4 + 0.2117342289834*a5 - 0.028417459543695*a6;
			mz1 = 0.2*mz2;
		}
		mz = (mz2*(mach - 0.0) + mz1 * (0.8 - mach)) / (0.8 - 0.0);

	}
	else {
		if (mach < 1.2) {
			if (alpha < 0) {
				mz1 = 0.30779142651144 - 0.277854038854256*a1 - 6.65707437924101*a2 - 10.3248375924987*a3 - 6.78103257929276*a4 - 2.30964522348498*a5 - 0.401738062965339*a6 - 0.0279362656081643*a7;
				mz2 = 0.305132030299153 + 0.397635804317453*a1 - 3.28505883652454*a2 - 5.17673270676827*a3 - 3.3520987433414*a4 - 1.22911795118107*a5 - 0.252937913593531*a6 - 0.0219642252350793*a7;

			}
			else {
				mz1 = 0.344548567530313 + 0.562705635856*a1 - 0.851458363904556*a2 + 0.610883020728777*a3 - 0.5142261300719*a4 + 0.2117342289834*a5 - 0.028417459543695*a6;
				mz2 = 0.326143039578254 + 0.838396408719481*a1 - 1.55964865562381*a2 + 1.20078969161253*a3 - 0.622928144295143*a4 + 0.169562746862479*a5 - 0.0166681591427948*a6;

			}
			mz = (mz2*(mach - 0.8) + mz1 * (1.2 - mach)) / (1.2 - 0.8);
		}
		else {
			if (mach < 2.0) {
				if (alpha < 0) {
					mz1 = 0.305132030299153 + 0.397635804317453*a1 - 3.28505883652454*a2 - 5.17673270676827*a3 - 3.3520987433414*a4 - 1.22911795118107*a5 - 0.252937913593531*a6 - 0.0219642252350793*a7;
					mz2 = 0.133393433017781 + 0.189771162177409*a1 + 0.294054590699251*a2 + 3.8738660326121*a3 + 5.4457292732551*a4 + 2.93842597002652*a5 + 0.70726993592008*a6 + 0.0641323809569589*a7;
				}
				else {
					mz1 = 0.326143039578254 + 0.838396408719481*a1 - 1.55964865562381*a2 + 1.20078969161253*a3 - 0.622928144295143*a4 + 0.169562746862479*a5 - 0.0166681591427948*a6;
					mz2 = 0.146009887776836 + 0.452070060631612*a1 + 0.151952550632889*a2 - 0.83782786872091*a3 + 0.494552561082258*a4 - 0.116890177238211*a5 + 0.0108109220801909*a6;
				}
				mz = (mz2*(mach - 1.2) + mz1 * (2.0 - mach)) / (2.0 - 1.2);
				//mz = mz2;
			}
			else {
				if (mach < 4.5) {
					if (alpha < 0) {
						mz2 = 0.107346213765651 + 0.458250461406157*a1 + 1.55000415146429*a2 + 3.74337646269654*a3 + 3.46977770226172*a4 + 1.38534982386379*a5 + 0.235344911650548*a6 + 0.012506244693478*a7;

						mz1 = 0.133393433017781 + 0.189771162177409*a1 + 0.294054590699251*a2 + 3.8738660326121*a3 + 5.4457292732551*a4 + 2.93842597002652*a5 + 0.70726993592008*a6 + 0.0641323809569589*a7;
					}
					else {
						mz2 = 0.11842066129015 + 0.457582935321605*a1 - 0.273059225480615*a2 - 0.109033887276009*a3 - 0.0103977310363484*a4 + 0.0490413170058151*a5 - 0.00972577829860893*a6;
						mz1 = 0.146009887776836 + 0.452070060631612*a1 + 0.151952550632889*a2 - 0.83782786872091*a3 + 0.494552561082258*a4 - 0.116890177238211*a5 + 0.0108109220801909*a6;
					}
					mz = (mz2*(mach - 2.0) + mz1 * (4.5 - mach)) / (4.5 - 2.0);
				} else {
					if (alpha < 0) {
						mz = 0.107346213765651 + 0.458250461406157*a1 + 1.55000415146429*a2 + 3.74337646269654*a3 + 3.46977770226172*a4 + 1.38534982386379*a5 + 0.235344911650548*a6 + 0.012506244693478*a7;
					}
					else {
						mz = 0.11842066129015 + 0.457582935321605*a1 - 0.273059225480615*a2 - 0.109033887276009*a3 - 0.0103977310363484*a4 + 0.0490413170058151*a5 - 0.00972577829860893*a6;

					}
				}
			}
		}
	}
	return mz;
}




double MyOmegaY(double mach, double betta) {
	return 0;

}

double CxModel5(double mach, double alpha, double Hg)
{
	double cxNoMod = CxPas(mach, alpha, Hg);
	double cxNoMod2 = 0;
	double cxMod2 = 0;
	double a1 = alpha;
	double a2 = alpha * alpha;
	double a3 = a2 * a1;
	double a4 = a3 * a1;
	double a5 = a4 * a1;
	double a6 = a5 * a1;
		if (alpha < 0) {			
			cxNoMod2 = 0.365580708328933 + 0.597715976382403*a1 + 0.557140991020652*a2 + 2.52599088374278*a3 + 2.2059056136048*a4 + 0.729553688770536*a5 + 0.0855518761652196*a6;
			cxMod2 = -(0.867991965540824 + 2.06659449121246*a1 - 0.515078099190406*a2 - 0.685172699073381*a3 - 0.118204670577509*a4);
		}
		else {
			cxNoMod2 = 0.38349254408567 + 0.529713035341652*a1 + 1.2050290101389*a2 - 1.45954297097336*a3 + 0.466491292961852*a4 - 0.0486311373721188*a5;
			cxMod2 = -(0.834571510452623 + 0.764862047837983*a1 - 0.120483693237789*a2 - 0.313304863237796*a3 + 0.0690238935156829*a4);
		}
		double k = (abs(cxMod2) > 1E-3) ? cxNoMod / cxMod2 : 2;
		double cx = (abs(k) > 2) ? 2 * cxMod2 : k * cxMod2;
		//double cx = cxNoMod * cxMod2 / cxNoMod2;
		
	return cx;
}

double CyModel5(double mach, double alpha, double Hg)
{
	double cyNoMod = CyAlPas(mach, alpha, Hg);
	double cyNoMod2 = 0;
	double cyMod2 = 0;
	double a1 = alpha;
	double a2 = alpha * alpha;
	double a3 = a2 * a1;
	double a4 = a3 * a1;
	double a5 = a4 * a1;
	double a6 = a5 * a1;
	if (alpha < 0) {
		cyNoMod2 = 0.740036659580577 + 5.76397325556497*a1 - 8.23184438631578*a2 - 9.35855376619993*a3 - 3.20462594417264*a4 - 0.389311584226631*a5;
		cyMod2 = -(-0.726519155970733 - 10.5516308248263*a1 - 2.66837005568588*a2 + 0.18462187927201*a3 + 0.0132927379145718*a4);
	}
	else {
		cyNoMod2 = 0.441581999540912 + 2.46220116703884*a1 + 9.29525221217141*a2 - 10.70759758376*a3 + 4.55991812468966*a4 - 0.914541243489724*a5 + 0.0679474605594143*a6;
		cyMod2 = -(-0.601745550772202 - 7.7572812656523*a1 + 2.69965775130698*a2);
	}
	double k = (abs(cyMod2)>1E-3)? cyNoMod / cyMod2: 2;
	double cy = (abs(k)>2)? 2 * cyMod2 : k * cyMod2;
	
	return cy;
}

double MzModel5(double mach, double alpha, double Hg)
{
	double mzNoMod = MzAlphaPas(mach, alpha, Hg);
	double mzNoMod2 = 0;
	double mzMod2 = 0;
	double a1 = alpha;
	double a2 = alpha * alpha;
	double a3 = a2 * a1;
	double a4 = a3 * a1;
	double a5 = a4 * a1;
	double a6 = a5 * a1;
	double a7 = a6 * a1;
	if (alpha < 0) {
		mzNoMod2 = 0.133393433017781 + 0.189771162177409*a1 + 0.294054590699251*a2 + 3.8738660326121*a3 + 5.4457292732551*a4 + 2.93842597002652*a5 + 0.70726993592008*a6 + 0.0641323809569589*a7;
	}
	else {
		mzNoMod2 = 0.146009887776836 + 0.452070060631612*a1 + 0.151952550632889*a2 - 0.83782786872091*a3 + 0.494552561082258*a4 - 0.116890177238211*a5 + 0.0108109220801909*a6;
	}
	if (alpha*toDeg < -30) {
		mzMod2 = 0.558383519820216 - 0.583279301932279*a1 - 4.01963134701603*a2 - 3.87938841882927*a3 - 1.34627688517851*a4 - 0.156749554196478*a5;
	}
	else {
		mzMod2 = 0.404276892052999 + 0.166787858027403*a1 - 0.341444013885305*a2 - 0.0050931963029156*a3 + 0.0249608605231202*a4;
	}
	double k = (abs(mzMod2) > 1E-3) ? mzNoMod / mzMod2 : 1.5;
	if (k < 0) { k = -k; }
	double mz = (abs(k) > 1.5) ? 1.5 * mzMod2 : k * mzMod2;
	//double mz = mzNoMod * mzMod2 / mzNoMod2;
	
	return mz;
}

double CxModel8(double mach, double alpha, double Hg)
{
	double cxNoMod = CxPas(mach, alpha, Hg);
	double cxNoMod2 = 0;
	double cxMod2 = 0;
	double a1 = alpha;
	double a2 = alpha * alpha;
	double a3 = a2 * a1;
	double a4 = a3 * a1;
	double a5 = a4 * a1;
	double a6 = a5 * a1;
	if (alpha < 0) {
		cxNoMod2 = 0.365580708328933 + 0.597715976382403*a1 + 0.557140991020652*a2 + 2.52599088374278*a3 + 2.2059056136048*a4 + 0.729553688770536*a5 + 0.0855518761652196*a6;
		cxMod2 = 1.60083391041373 + 2.84327704113597*a1 - 0.745037275188489*a2 - 0.909576231555973*a3 - 0.157003027056968*a4;
	}
	else {
		cxNoMod2 = 0.38349254408567 + 0.529713035341652*a1 + 1.2050290101389*a2 - 1.45954297097336*a3 + 0.466491292961852*a4 - 0.0486311373721188*a5;
		cxMod2 = 1.55718119489783 + 0.929250741679049*a1 - 0.2507148706551*a2 - 0.52608418507473*a3 + 0.129129984418612*a4;
	}
	double k = (abs(cxMod2) > 1E-3) ? cxNoMod / cxMod2 : 2;
	double cx = (abs(k) > 2) ? 2 * cxMod2 : k * cxMod2;
	//double cx = cxNoMod * cxMod2 / cxNoMod2;

	return cx;
}

double CyModel8(double mach, double alpha, double Hg)
{
	double cyNoMod = CyAlPas(mach, alpha, Hg);
	double cyNoMod2 = 0;
	double cyMod2 = 0;
	double a1 = alpha;
	double a2 = alpha * alpha;
	double a3 = a2 * a1;
	double a4 = a3 * a1;
	double a5 = a4 * a1;
	double a6 = a5 * a1;
	if (alpha < 0) {
		cyNoMod2 = 0.740036659580577 + 5.76397325556497*a1 - 8.23184438631578*a2 - 9.35855376619993*a3 - 3.20462594417264*a4 - 0.389311584226631*a5;
	//	cyMod2 = -(-0.726519155970733 - 10.5516308248263*a1 - 2.66837005568588*a2 + 0.18462187927201*a3 + 0.0132927379145718*a4);
	}
	else {
		cyNoMod2 = 0.441581999540912 + 2.46220116703884*a1 + 9.29525221217141*a2 - 10.70759758376*a3 + 4.55991812468966*a4 - 0.914541243489724*a5 + 0.0679474605594143*a6;
	//	cyMod2 = -(-0.601745550772202 - 7.7572812656523*a1 + 2.69965775130698*a2);
	}
	cyMod2 = -(-0.467668354611325 - 7.33812913442053*a1 + 0.678596882716044*a2 + 1.23824644824045*a3 - 0.0468748000871884*a4 - 0.0502748312806952*a5);
	double k = (abs(cyMod2) > 1E-3) ? cyNoMod / cyMod2 : 2;
	double cy = (abs(k) > 2) ? 2 * cyMod2 : k * cyMod2;

	return cy;
}

double MzModel8(double mach, double alpha, double Hg)
{
	double mzNoMod = MzAlphaPas(mach, alpha, Hg);
	double mzNoMod2 = 0;
	double mzMod2 = 0;
	double a1 = alpha;
	double a2 = alpha * alpha;
	double a3 = a2 * a1;
	double a4 = a3 * a1;
	double a5 = a4 * a1;
	double a6 = a5 * a1;
	double a7 = a6 * a1;
	if (alpha < 0) {
		mzNoMod2 = 0.133393433017781 + 0.189771162177409*a1 + 0.294054590699251*a2 + 3.8738660326121*a3 + 5.4457292732551*a4 + 2.93842597002652*a5 + 0.70726993592008*a6 + 0.0641323809569589*a7;
	}
	else {
		mzNoMod2 = 0.146009887776836 + 0.452070060631612*a1 + 0.151952550632889*a2 - 0.83782786872091*a3 + 0.494552561082258*a4 - 0.116890177238211*a5 + 0.0108109220801909*a6;
	}
	if (alpha*toDeg < -30) {
		mzMod2 = 0.891858901447927 - 0.132019770510859*a1 - 3.90302811806877*a2 - 4.03910285082207*a3 - 1.46895053774131*a4 - 0.179368830244217*a5;
	}
	else {
		mzMod2 = 0.51617679323505 + 0.043286804339664*a1 - 0.438616569182197*a2 + 0.112336744223521*a3;
	}
	double k = (abs(mzMod2) > 1E-3) ? mzNoMod / mzMod2 : 1.5;
	if (k < 0) { k = -k; }
	double mz = (abs(k) > 1.5) ? 1.5 * mzMod2 : k * mzMod2;
	//double mz = mzNoMod * mzMod2 / mzNoMod2;

	return mz;
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
			if (alphaDeg < 181.25) {
				my1 = 1.177265E-3 + 2.0660087E-3*alphaDeg + 5.1360683E-5*alphaDeg2 - 1.99451668E-6*alphaDeg3 + 1.61694008413613E-8*alphaDeg4 - 3.90498782894E-11*alphaDeg5;
			}
			else {				
					my1 = -6.59739187 + 0.13780116*alphaDeg - 1.151492833125E-3*alphaDeg2 + 4.801379052E-6*alphaDeg3 - 9.9470177003001E-9*alphaDeg4 + 8.15104025386042E-12*alphaDeg5;				
			}
			if (alphaDeg < 161.818) {
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
				if (alphaDeg < 161.818) {
					my1 = 1.293096395E-3 + 3.4148346676E-3*alphaDeg + 4.7291469E-5*alphaDeg2 - 2.4200543938E-6*alphaDeg3 + 2.017695332E-8*alphaDeg4 - 4.85801728813182E-11*alphaDeg5;
				}
				else {
					
						my1 = 6.0075128 - 0.1103382468*alphaDeg + 7.33545E-4*alphaDeg2 - 2.089793715E-6*alphaDeg3 + 2.171407611E-9*alphaDeg4 - 5.2208322772E-14*alphaDeg5;					
				}


				if (alphaDeg < 180.625) {
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
	double alphaDeg5 = alphaDeg2 * alphaDeg3;
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
			if (alphaDeg < 180) {
				mx1 = -4.11595338E-3 + 4.322027523E-4*alphaDeg - 1.518620114E-4*alphaDeg2 + 2.3259526793E-6*alphaDeg3 - 1.2024717703E-8*alphaDeg4 + 2.06506966141321E-11*alphaDeg5;
			}
			else {
				mx1 = 6.873320132 - 0.1249343665*alphaDeg + 8.830828935E-4*alphaDeg2 - 3.04725145156E-6*alphaDeg3 + 5.17817218653E-9*alphaDeg4 - 3.49735022619433E-12*alphaDeg5;
			}
			if (alphaDeg < 178.812) {
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
				if (alphaDeg < 178.812) {
					mx1 = -4.819184741E-3 + 3.720456549E-4*alphaDeg - 2.3052581254E-4*alphaDeg2 + 3.29557444863E-6*alphaDeg3 - 1.529297242E-8*alphaDeg4 + 2.24075890939032E-11*alphaDeg5;
				}
				else {

					mx1 = 13.56323183 - 0.22266758673*alphaDeg + 1.34401620049E-3*alphaDeg2 - 3.61050893514E-6*alphaDeg3 + 4.0458625237E-9*alphaDeg4 - 1.1729417603405E-12*alphaDeg5;

				}


				if (alphaDeg < 180) {
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
