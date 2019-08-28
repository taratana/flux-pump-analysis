/*
	This program is for analysis of flux pump using pancake coil.
	In addition, self of mutual inductance, contact resistance, etc. are all reasonably calculated.
	This program is made by Mato in 2019.
*/

#define _CRT_SECURE_NO_WARNINGS
#include"Circuit.h"
#include "Inductance.h"
using namespace Eigen;
using namespace std;

class StopWatch {
private:
	time_t  st_t, fin_t;  /*時間計測用（実時間） */
	clock_t st_ct, fin_ct; /*時間計測用（CPU時間）*/

public:
	void start() { /*時間計測開始*/
		st_ct = clock();
		time(&st_t);
	}
	void stop() { /*時間計測終了*/
		fin_ct = clock();
		time(&fin_t);
	}
	double getTime() {
		return difftime(fin_t, st_t);
	}
};


int main() {
	StopWatch record_time;
	record_time.start();


	//Inductance inductance;
	//inductance.sc1->setCoilParameter(0.0875, 0.1125, -0.0425, -0.0175, 200);
	//inductance.sc2->setCoilParameter(0.0875, 0.1125, 0.0175, 0.0425, 200);
	//cout << inductance.M() << endl;
	//MatrixXd M_mat(NOP, NOP);
	//inductance.CalcInductance(M_mat);

	//Circuit *circuit;
	//circuit.CalcM_ij();
	//circuit.CalcC(0, 1);

	//for (int i = 1; i <= 1; i++) {
	//	circuit = new Circuit();
	//	circuit->setNOP(i);
	//	circuit->CalcCircuit("NOP"+to_string(i));
	//	delete circuit;
	//}


	static const int N = 200+1;
	Circuit *circuit[N];
	char filename[100];
	double conv_info[N][3];
	double filter_initial = 0.8e-6;//  0.0004955;
	double filter_final = 3.0e-1;
	double log_space_factor = pow(filter_final / filter_initial, 1.0 / (N - 1));
	for(int j=12;j<=12;j++){
#pragma omp parallel for schedule (guided,1)
		for (int i = 0; i < N; i++) {
			sprintf(filename, "ExpData/09/%02d/Data%d", j, i);
			circuit[i] = new Circuit();
			circuit[i]->setNOP(j);
			if (i == 0) {
				conv_info[i][0] = 0;
			} else {
				conv_info[i][0] = filter_initial*pow(log_space_factor, i-1);

			}
			circuit[i]->setFilterInductance(conv_info[i][0]);
			circuit[i]->CalcCircuit(filename);
			circuit[i]->getConvergenceInfo(&conv_info[i][1], &conv_info[i][2]);
			delete circuit[i];
			cout << "Circuit" << i << " Finished." << endl;
		}
		sprintf(filename, "ExpData/09/%02d/conv-Lf.csv", j);
		ofstream output_Lf(filename);
		for (int i = 0; i < N; i++)
			output_Lf << conv_info[i][0] << "," << conv_info[i][1] << "," << conv_info[i][2] << endl;
		output_Lf.close();
	}


//	static const int N = 1;
//	Circuit *circuit[N];

//	double conv_info[N][3];
//#pragma omp parallel for schedule (dynamic,1)
//	for (int i = 0; i < N; i++) {
//		circuit[i] = new Circuit();
//		conv_info[i][0] = 0;
//		//circuit[i]->setR_dyn(conv_info[i][0]);
//		circuit[i]->setFilterInductance(conv_info[i][0]);
//		circuit[i]->CalcCircuitWithoutfilter("ExpData/01/Data" + to_string(i + 1));
//		circuit[i]->getConvergenceInfo(&conv_info[i][1], &conv_info[i][2]);
//		delete circuit[i];
//		cout << "Circuit" << i << " Finished." << endl;
//	}
//	ofstream output_Lf("ExpData/01/conv-Lf.csv");
//	for (int i = 0; i < N; i++)
//		output_Lf << conv_info[i][0] << "," << conv_info[i][1] << "," << conv_info[i][2] << endl;
//	output_Lf.close();


	record_time.stop();
	double elapsed = record_time.getTime();
	cout << "Time Elapsed:" << endl << (int)elapsed / 3600 << "hrs " << ((int)elapsed % 3600) / 60.0 << "min" << endl;

}