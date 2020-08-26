#include "pch.h"
#include "Atmosphere.h"

//�����������
Atmosphere::Atmosphere() {
	fillArrPandT();
	i = 0;
}




// ������� �� �������� ���������������� ������ ���������� ���� = dT/dH
double Atmosphere::betaKoef(int i) {
	
	return (Beta[1][i]);
}
int Atmosphere::index(double geoPotH) {
	//�������� ���������������� ������ - ������ ������� ������� � �
	//�������� ������������ ����=dT/dH - ������ ������ �������
	//double gH[12] = { -2000.0, 0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 85000.0, 94.0E+3, 102.45E+3, 117.777E+3 };
	int i = 0;
	while ((geoPotH >= gPotHStand[i + 1])&(i < 11)) i++;

	return (i);
}
double Atmosphere::tFunc(double h) {  // h = �������������� ������
	double gh = gHFh(h);  //���������������� ������
	if ((gh < gPotHStand[i]) || (gh > gPotHStand[i + 1])) {
		i = index(gh);
	}
	
	double b = betaKoef(i);
	double t = temperature[i] + b * (gh - gPotHStand[i]);
	return t;
}


//���������
double Atmosphere::roFunc(double h) {
	double t = tFunc(h);
	double p = pFunc(h);
	double ro = p / (R*t);
	return (ro);
}
//�������� �����
double Atmosphere::aFunc(double h) {
	double t = tFunc(h);
	double a = 20.046796*sqrt(t);
	return (a);
}
// ������������ �������� �������
double Atmosphere::muFunc(double h) {
	double t = tFunc(h);
	double mu = BS * pow(t, 1.5) / (t + S);
	return (mu);
}
// �������������� ��������
double Atmosphere::nuFunc(double h) {
	double mu = muFunc(h);
	double ro = roFunc(h);
	double nu = mu / ro;
	return (nu);
}
// �������������� ������
double Atmosphere::hFromgH(double gH) {
	const double RZU = 6356767.0;
	double h = RZU * gH / (RZU - gH);
	return (h);
}
//���������������� ������
double Atmosphere::gHFh(double h) {
	const double RZU = 6356767.0;
	double gH = RZU * h / (RZU + h);
	return (gH);
}

void Atmosphere::fillArrPandT() {
	double b;
	for (int i = 0; i <= 12; i++) {

		switch (i) {
		case 0:

			b = betaKoef(0);
			temperature[0] = TST + b * (gPotHStand[0] - gPotHStand[1]);
			if (b == 0)
				pGostStand[0] = PST * exp(GST*(gPotHStand[1] - gPotHStand[0]) / (R * TST));
			else pGostStand[0] = PST * pow((1.0 + b * ((gPotHStand[1] - gPotHStand[0]) / temperature[0])), (GST / (b *R)));
			break;
		case 1:
			pGostStand[1] = PST;
			temperature[1] = TST;
			break;
		default:
			b = betaKoef(i - 1);
			temperature[i] = temperature[i - 1] + b * (gPotHStand[i] - gPotHStand[i - 1]);

			if (b == 0)
				pGostStand[i] = pGostStand[i - 1] * exp((-1)*GST*(gPotHStand[i] - gPotHStand[i - 1]) / (R*temperature[i]));
			//pGostStand[i] = pGostStand[i - 1] * pow(2.71829, ((-1)*GST*(gPotHStand[i] - gPotHStand[i - 1]) / (R*temperature[i])));
			else pGostStand[i] = pGostStand[i - 1] * pow((1 + b * ((gPotHStand[i] - gPotHStand[i - 1]) / temperature[i - 1])), ((-1)*GST / (b * R)));
		}
		//pGostStand[i] = roundBan1(pGostStand[i]);
	}
	return;
}
double Atmosphere::pFunc(double h) {
	double gh = gHFh(h);  //���������������� ������
	if ((gh < gPotHStand[i]) || (gh > gPotHStand[i + 1])) {
		i = index(gh);
	}
	double b = betaKoef(i);
	double tPr = temperature[i] + b * (gh - gPotHStand[i]);
	double pPr;
	if (b == 0)
		pPr = pGostStand[i] * exp((-1)*GST*(gh - gPotHStand[i]) / (R*tPr));
	else pPr = pGostStand[i] * pow((1 + b * ((gh - gPotHStand[i]) / temperature[i])), ((-1)*GST / (b * R)));
	return (pPr);
}

Atmosphere::~Atmosphere() {}



