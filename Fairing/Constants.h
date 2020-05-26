#pragma once
const double PI = 3.14159265358979323846264338328;
const double toDeg = 180 / PI;
const double toRad = PI / 180;

const int N = 6;

const double X0 = 0;
const double Y0 = 10000 - 150 * N;
const double Z0 = 0;

const double ALPHA0 = 0;
const double M0 = 1400;

const double W_X0 = 0;
const double W_Y0 = 0;
const double W_Z0 = 0;

const double I_X0 = 59110.61 * M0 / 40816.77;
const double I_Y0 = 327307.83 * M0 / 40816.77;
const double I_Z0 = 364828.28 * M0 / 40816.77;
const double I_XY0= 16389 * M0 / 40816.77;
const double I_YX0 = 16389 * M0 / 40816.77;

const double D_M = 0.95;
const double L = 7;

//const double V0 = 1200; //модуль скорости


const double H = 0.000025;
const double T_FIN = 10.0;
const double EPS1 = 1E-6;

const double KSI_SST = 0.35;
const double KSI_SSN = 0.35;
const double KSI_SSE = 0.35;

const double K_SST = 0.95;
const double K_SSN = 0.95;
const double T_SSE = 0.01;
/////////////////////////////////////////////////////////
const double OMEGA_EARTH = 7.2921158553E-5;

const double RA_EL = 6378137;
const double ALPHA_EL = 1 / 298.25784;

const double LAT = 45.966111*toRad;
const double LON = 63.307778*toRad;
const double AZIM = 61.3*toRad;



const double PITCH0 = (20) * toRad; //тангаж в радианах
const double ROLL0 = (0)*toRad; //угол крена в радианах
const double YAW0 = (0)*toRad; //угол рысканья в радианах

//const double L1 = 21.87;
//const double L_CONE = 4.4;
//const double D_m = 2.1 * 2;
//
//
//const double L_FULL = 58;
//const double XD_FULL = 35.04;
