/*
[How to use this program]
This program calculate self or mutual inductance.
Calculation only needs the shape of coils,
which is inner radius, outer radius, top coordinate, bottom coordinate.
Below is example of code:
	Inductance inductance;
	inductance.sc1->setCoilParameter(0.0875, 0.1125, -0.0425, -0.0175, 200);
	inductance.sc2->setCoilParameter(0.0875, 0.1125, 0.0175, 0.0425, 200);
	cout << inductance.M() << endl;
Then, you get mutual indactance.
If you want to calculate self indactance, you have to set same shape into sc1 and sc2.

For more detail about calculation, plese refer to this paper "Usable ranges of some expressions 
for calculation of the self-inductance of a circular coil of rectangular cross section."
This program is made by Kurauchi and modified by Mato in 2019.
*/


#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include "GaussLegendreQuadrature.h"
#include "headings.h"
#include<Eigen/Dense>


static const int N0 = 500;
static const int NR = 500;
static const int NZ = 500;
static const double REBCO_WIDTH = 4e-3;
static const double INSULATOR_WIDTH = 1e-3;
static const int GAUSS_POINTS_NUMBER = 20;


/*
矩形コイルの周りの磁界計算用
*/
class Inductance {
private:
	class SingleCoil {
	public:
		double a1, a2, z1, z2;
		double J;
		double Rp, Zp;
		double n0, nr, nz;
		int number_of_turn;
		Rosetta::GaussLegendreQuadrature<GAUSS_POINTS_NUMBER> gaussIntegration;

		SingleCoil();

		/*磁場を計算する*/
		double Bz(double Rp, double Zp);
		double Hr();
		double Hz();

		double Hz0(double a, double z);
		double Hr1(double a, double z);
		double Hr2(double a, double z);
		double Hr2_1(double a, double z);
		double Hr2_2(double a, double z);
		double Hr2_3(double a, double z);
		double Hz1(double a, double z);
		double Hz1_1(double a, double z);
		double Hz1_2(double a, double z);
		double Hz2(double a, double z);
		double Hz2_1(double a, double z);
		double Hz2_2(double a, double z);
		double Hz2_3(double a, double z);
		double Hz3(double a, double z);
		double R(double a, double z, double theta);
		double CalcJ();
		void setRpZp(double _Rp, double _Zp);
		void setDivNum(int n0, int nr, int nz);
		void setCoilParameter(double _a1, double _a2, double _z1, double _z2, int turn);
	};


public:
	//sc1 for the coil generating magnetic field 
	SingleCoil *sc1, *sc2;

	Inductance();
	void CalcInductance(Eigen::MatrixXd &M);
	double M();
	double ri(int i);
	double rj(int j);
	double zh(int h);

};


