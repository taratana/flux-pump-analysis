#pragma once
#include "headings.h"
#include "DataFile.h"
#include "Inductance.h"
#include "GMRES.h"
using namespace Eigen;



class Circuit {
private:
	double i_2_freq, i_2_mag;
	double L_f;
	double total_inductance;
	double b_app_rate;
	double b_mag, b_freq;
	double b_app_length, b_app_width;
	double r_dyn_c;
	double r_dyn_mag;
	int nop;
	int nou;
	double convergence_time;
	double convergence_value;

	Inductance inductance;

	VectorXd time;
	VectorXd I_L, I_B;
	VectorXd C;		//C for KI=C
	VectorXd R_c;	//Contact resistance
	MatrixXd I_theta;	//NOPxTIME_DIV_NUM
	MatrixXd I_r;		//NOPxTIME_DIV_NUM
	MatrixXd K;		//Coefficient of circuit equation in the form of matrix
	MatrixXd M_ij;	//Self or mutual inductance

	class DataFile *opd;
public:
	Circuit();
	void InitMatrices();
	void setI2(double i_2_mag, double i_2_freq);
	void setp(double p);
	void setFilterInductance(double L_f);
	void setR_dyn(double r_dyn);
	void setR_c(double r_c);
	void setResistance(double R_c, double R_joint);
	void setNOP(int);
	void getConvergenceInfo(double *conv_t, double *conv_v);

	void CalcM_ij();
	void CalcK(double t);
	void CalcC(double t, int i_for_time);
	void CalcRc();
	double CalcR_dyn(double t);
	double CalcI_2(double t);
	void CalcCircuit(string filename);
	void CalcConvergenceInfo();
	void PrintInfo(string filename);
};