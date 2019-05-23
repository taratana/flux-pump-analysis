#pragma once
#include "constants.h"
#include<math.h>
#include<fstream>
#include<omp.h>
#include<iostream>
#include<iomanip>
#include<string>

class circuit {
private:
	double r_c, r_joint;
	double i_2_freq, i_2_mag;
	double l_coil, l_filter;
	double b_app_rate;
	double b_mag, b_freq;
	double b_app_length, b_app_width;
	double r_dyn_c;
	double r_dyn_mag;
public:
	circuit();
	void setI2(double i_2_mag,double i_2_freq);
	void setp(double p);
	void setRdyn(double rdyn);
	void CalcNoFilterCircuit(std::string output_filename);
	void CalcFilterExistCircuit(std::string output_filename);
	double CalcR_dyn(double t);
	double CalcI_2(double t);
	double FuncForNoFilterCircuit(double t, double i_theta);
	double FuncForFilterCircuit1(double t, double i_theta, double i_theta_d);
	double FuncForFilterCircuit2(double t, double i_theta, double i_theta_d);

};