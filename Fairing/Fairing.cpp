// Fairing.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include "Runge.h"
#include "Launch.h"
#include "Drop.h"
#include "DistrFull.h"
	int flag;
	void dropInteg(Drop &Rocket2, Drop &buf) {
		double h1 = H;
		int flag1 = 0;
		while ((abs(Rocket2.getH()) > 1E-1) && (Rocket2.getParam().vect[13] < 1000)) {
			if (Rocket2.getH() < 1E+3) {
				buf = Rocket2;
			}
			if (flag1 == 0) {
				h1 = H / pow(2, Rocket2.getParam().vect[8]);
				if (h1 < 1E-8) {
					h1 = 1E-8;
				}
			}
			//h1 = H / pow(2, Rocket2.getParam().vect[8]);
			//	traj2.push_back(Rocket2);
			runge(Rocket2, h1);
			if (Rocket2.getH() < 0) {
				Rocket2 = buf;
				h1 /= 2;
				flag1 = 1;
			}
		}
		if (Rocket2.getParam().vect[13] > 1000) {
			flag = 1;
		}
	}

	void dropInteg(Drop &Rocket2, Drop &buf, std::ofstream &fout) {
		double h1 = H;
		int count = 0;
		while ((abs(Rocket2.getH()) > 1E-1) && (Rocket2.getParam().vect[13] < 5000)) {
			if (Rocket2.getH() < 1E+3) {
				buf = Rocket2;
			}
			if (count % 100 == 0/*(fabs(Rocket3.getParam().vect[13] * 1000 - round(Rocket3.getParam().vect[13] * 1000)) < EPS1)
					&& ((int)round(Rocket3.getParam().vect[13] * 1000) % 5 == 0)*/) {
				Rocket2.printParam(fout);

			}
			//	traj2.push_back(Rocket2);
			runge(Rocket2, h1);
			count++;
		//	Rocket2.printParam(fout);
			if (Rocket2.getH() < 0) {
				Rocket2 = buf;

				h1 /= 2;
			}
		}
		Rocket2.printParam(fout);
		if (Rocket2.getParam().vect[13] > 1000) {
			flag = 1;
		}
		
	}

	int main()
	{

		
	
		double r = RA_EL * (1 - ALPHA_EL * sin(LAT)*sin(LAT));
		double Vel = 2000;
		double V1 = Vel * cos(PITCH0) * cos(YAW0);
		double V2 = Vel * sin(PITCH0);
		double V3 = Vel * cos(PITCH0) * sin(YAW0);
		std::vector <double> b = { V1,V2,V3,/**/0,r + 80E+3,0,/**/0,0,-1,/**/PITCH0 , YAW0, ROLL0,0,M0 };
		
		Drop buf3(b);
		Drop Rocket3(b);

		//ПОМЕНЯТЬ НА ТОЧНЫЕ ЗНАЧЕНИЯ
		double Ix1 = I_X0;
		double Iy1 = I_Y0;
		double Iz1 = I_Z0;
		double Smm = 10.746/2;
		double Ll = 9.38;
			Rocket3.setStageParam(30000, 0, 0, Smm, 0, Ll, 0, Ix1, Iy1, Iz1);
			Rocket3.addErotation(LAT, AZIM);
			Rocket3.nonIntegr();
			buf3.setStageParam(30000, 0, 0, Smm, 0, Ll, 0, Ix1, Iy1, Iz1);
			buf3.addErotation(LAT, AZIM);
			buf3.nonIntegr();

			std::string bfilename = "varStepRunge.txt";
			std::ofstream fout2(bfilename);


			double h1 = H;

			int count = 0;
			while (Rocket3.getH() > 0) {
				h1 = H / pow(2, Rocket3.getParam().vect[8]);
				if (count % 800 == 0/*(fabs(Rocket3.getParam().vect[13] * 1000 - round(Rocket3.getParam().vect[13] * 1000)) < EPS1)
					&& ((int)round(Rocket3.getParam().vect[13] * 1000) % 5 == 0)*/) {
					Rocket3.printParam(fout2);

				}
				count++;
			//	traj.push_back(Rocket3);
				if (Rocket3.getH() < 1E+3) {
					buf3 = Rocket3;
					//buf = Rocket;
				}
			
				if (count == 31600) {
					double g = 0;
				}
				runge(Rocket3, h1);

			}
			Rocket3 = buf3;
			//Rocket = buf;
			//h1 = H;

			while (abs(Rocket3.getH() > 1E-1)) {
				buf3 = Rocket3;
				runge(Rocket3, h1);
				if (Rocket3.getH() < 0) {
					Rocket3 = buf3;
					h1 /= 2;
				}
			}
			Rocket3.printParam(fout2);

			/*std::string bfilename1 = "pr.txt";
			std::ofstream fout3(bfilename1);
*/

			std::string bfilename3 = "d_paramVarStepVarAtm.txt";
			std::ofstream fout3(bfilename3);

			std::string bfilename4 = "res_dropVarStepVarAtm.txt";
			std::ofstream fout4(bfilename4);
			fout4 << '\t';
			Rocket3.printParam(fout4);
			std::random_device rd;
			std::mt19937 gen(rd());
			std::normal_distribution<> dis(0, 1.0*0.3);

			double dV = Vel * 0.05;
			double dWz = 1 * 0.05;
			double dPitch = PITCH0 * 0.05;
			double rund = dis(gen);
		//	#pragma loop(hint_parallel(8));
			for (int i = 0; i < 100; i++) {
				//std::string count = std::to_string(i);
				//std::string filename = "res_varAtm" + std::to_string(1)+".txt";
				//std::ofstream fout24(filename);
				rund = dis(gen);
				double dVz = /*dV * 0.5 * dis(gen)*/0;
				double dVy = /*sqrt(dV*dV - dVz * dVz) * dis(gen)*/0;
				double dVx = /*sqrt(dV*dV - dVz * dVz - dVy * dVy)*dis(gen)*/0;
				double dWwz = /*dWz * dis(gen)*/0;
				double dWwy = 0/*1E-1 * dis(gen)*/;
				double dWwx = 0/*1E-1 * dis(gen)*/;
				double ddPitch = 0/*dPitch * dis(gen)*/;
				double windX = 0 * dis(gen);
				double windZ = sqrt(0 * 0 - windX * windX)/**dis(gen)*/;


				b[9] += ddPitch;
				b[2] += dVz;
				b[1] += dVy;
				b[0] += dVx;
				b[8] += dWwz;
				b[6] += dWwx;
				b[7] += dWwy;

				/*Drop buf3(b);
				Drop Rocket3(b);*/
				distrFull buf3(b, windX, windZ);
				distrFull Rocket3(b, windX, windZ);

				Rocket3.setStageParam(30000, 0, 0, Smm, 0, Ll, 0, Ix1, Iy1, Iz1);
				Rocket3.addErotation(LON, AZIM);

				Rocket3.nonIntegr();
				buf3.setStageParam(30000, 0, 0, Smm, 0, Ll, 0, Ix1, Iy1, Iz1);
				buf3.addErotation(LON, AZIM);
				buf3.nonIntegr();



				dropInteg(Rocket3, buf3);
				//static int i = 0;
				//i++;
				std::cout << i;
				fout3 << i << '\t' << dVx << '\t' << dVy << '\t' << dVz << '\t' << ddPitch << '\t' << dWwz << '\t' << windX << '\t' << windZ << '\t' << '\t' << b[0] << '\t' << b[1] << '\t' << b[2] << '\t' << ddPitch << '\t' << dWwz << '\t' << windX << '\t' << windZ << '\n';
				fout4 << i << '\t';
				if (flag != 1) {
					Rocket3.printParam(fout4);
				}
				else {
					fout4 << i << '\n';
					flag = 0;
				}
				b[9] -= ddPitch;
				b[2] -= dVz;
				b[1] -= dVy;
				b[0] -= dVx;
				b[8] -= dWwz;

			
			}
		//for (int i = 0; i < 100; i++) {
		//	fout3 << 4.5*i/100 << '\t';
		//	fout3 << CxPas(4.5*i/100, 60*toRad, 5) << '\t';
		//	/*fout3 << CxPas(1.2, i*toRad, 5) << '\t';
		//	fout3 << CxPas(2.0, i*toRad, 5) << '\t';
		//	fout3 << CxPas(4.5, i*toRad, 5) << '\t' << '\t';*/

		//	fout3 << CyAlPas(4.5*i / 100, 60 * toRad, 5) << '\t';
		//	/*fout3 << CyAlPas(1.2, i*toRad, 5) << '\t';
		//	fout3 << CyAlPas(2.0, i*toRad, 5) << '\t';
		//	fout3 << CyAlPas(4.5, i*toRad, 5) << '\t' << '\t';*/

		//	//fout3 << MzAlphaPas(0.8, i*toRad, 5) << '\t';
		//	//fout3 << MzAlphaPas(1.2, i*toRad, 5) << '\t';
		//	//fout3 << MzAlphaPas(2.0, i*toRad, 5) << '\t';
		//	fout3 << MzAlphaPas(4.5*i / 100, 60 * toRad, 5) << '\t' << '\n';
		//}

		 std::cout << "Hello World!\n"; 

	}

