#pragma once
#include "headings.h"
#include "OutputData.h"

using namespace Eigen;

class Circuit {
private:
	double i_2_freq, i_2_mag;
	double L_f;
	double b_app_rate;
	double b_mag, b_freq;
	double b_app_length, b_app_width;
	double r_dyn_c;
	double r_dyn_mag;
	

	VectorXd time;
	VectorXd I_L, I_B;
	VectorXd C;		//C for KI=C
	VectorXd R_c;	//Contact resistance
	MatrixXd I_theta;	//NOPxTIME_DIV_NUM
	MatrixXd I_r;		//NOPxTIME_DIV_NUM
	MatrixXd K;		//Coefficient of circuit equation in the form of matrix
	MatrixXd M_ij;	//Self or mutual inductance

	class OutputData *opd;
	
public:
	Circuit();
	void Init();
	void setI2(double i_2_mag,double i_2_freq);
	void setp(double p);
	void setInductance(double L_c,double L_f);
	void setR_dyn(double r_dyn);
	void setResistance(double R_c,double R_joint);

	void CalcM_ij();
	void CalcK(double t);
	void CalcC(double t,int i_for_time);
	double CalcR_dyn(double t);
	double CalcI_2(double t);
	void CalcCircuit();



};