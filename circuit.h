#pragma once
#include "constants.h"
#include<math.h>
#include<fstream>
#include<omp.h>
#include<iostream>
#include<iomanip>

class circuit {
private:

public:
	void CalcNoFilterCircuit();
	void CalcFilterExistCircuit();
	double CalcR_dyn(double t);
	double CalcI_2(double t);
	double FuncForNoFilterCircuit(double t, double i_theta);
	double FuncForFilterCircuit1(double t,double i_theta, double i_theta_d);
	double FuncForFilterCircuit2(double t, double i_theta,double i_theta_d);

};