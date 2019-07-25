/*
	This program is for analysis of flux pump using pancake coil.
	In addition, self of mutual inductance, contact resistance, etc. are all reasonably calculated.
	This program is made by Mato in 2019.
*/

#define _CRT_SECURE_NO_WARNINGS
#pragma once
#include"circuit.h"
#include <Eigen/Dense>
#include "GMRES.h"
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


	//int n = 5;
	//VectorXd x = VectorXd::Zero(n);
	//x << 100,100,100,100,100;
	//VectorXd b(n);
	//b << -53,1,29,22,-5;

	//MatrixXd A(n, n);
	//A << 1, -2, 3, 4, -5,
	//	1, 1, 1, 1, 1,
	//	1, -2, 3, -4, 5,
	//	2, -1, -3, -2, 1,
	//	-1, -1, 1, 1, 1;


	//x = GMRES(A, b, x, 1e-9, 100, 10);

	//for (int i = 0; i < n; i++) {
	//	cout << x(i) << endl;
	//}
	//while (1);

	Inductance inductance;
	inductance.sc1->setCoilParameter(0.0875, 0.1125, -0.0425, -0.0175, 200);
	inductance.sc2->setCoilParameter(0.0875, 0.1125, 0.0175, 0.0425, 200);
	cout << inductance.M()*1000 << endl;
	


	record_time.stop();
	double elapsed = record_time.getTime();
	cout << "Time Elapsed:" << endl << (int)elapsed / 3600 << "hrs " << ((int)elapsed % 3600) / 60.0 << "min" << endl;
	
}