#include "circuit.h"

using namespace std;
using namespace Eigen;

Circuit::Circuit() {
	i_2_freq = i_2_FREQUENCY;
	i_2_mag = i_2_MAGNITUDE;
	L_f = FILTERING_COIL_INDUCTANCE;
	b_app_rate = B_APP_RATE;
	b_mag = B_MAGNITUDE;
	b_freq = B_FRECUENCY;
	b_app_length = B_APP_LENGTH;
	b_app_width = B_APP_WIDTH;
	r_dyn_c = C_OF_DYNAMIC_RESISTENCE;
	r_dyn_mag = 2 * b_app_length*b_app_width*b_freq*(b_mag + r_dyn_c * b_mag*b_mag) / CRITICAL_CURRENT_OF_YBCO_TAPE;

	Init();

}

//Initializing arrays
void Circuit::Init() {
	time = VectorXd::Zero(TIME_DIV_NUM);
	I_L = VectorXd::Zero(TIME_DIV_NUM);
	I_B = VectorXd::Zero(TIME_DIV_NUM);
	R_c = VectorXd::Zero(TIME_DIV_NUM);
	I_theta = MatrixXd::Zero(NOP, TIME_DIV_NUM);
	I_r = MatrixXd::Zero(NOP, TIME_DIV_NUM);

	C = VectorXd::Zero(NUMBER_OF_UNKOWNS, NUMBER_OF_UNKOWNS);
	K = MatrixXd::Zero(NUMBER_OF_UNKOWNS, NUMBER_OF_UNKOWNS);
	M_ij = MatrixXd::Zero(NOP, NOP);
}

void Circuit::setI2(double i_2_mag, double i_2_freq) {
	this->i_2_mag = i_2_mag;
	this->i_2_freq = i_2_freq;
}

void Circuit::setp(double p) {
	this->b_app_rate = p;
}

void Circuit::setInductance(double L_c, double L_f) {
	this->L_f = L_f;
}

void Circuit::setR_dyn(double r_dyn) {
	this->r_dyn_mag = r_dyn;
}

void Circuit::setResistance(double R_c, double R_joint) {

}

void Circuit::CalcM_ij() {
}

void Circuit::CalcK(double t) {
	double dt = TIME_TICK;

	K.setZero(NUMBER_OF_UNKOWNS, NUMBER_OF_UNKOWNS);

	K.block(0, 0, NOP, NOP) = M_ij;
	K.block(NOP, 0, NOP, NOP) = MatrixXd::Identity(NOP, NOP);
	K.block(0, NOP, NOP, NOP) = -1*R_c.asDiagonal();
	K.block(NOP, NOP, NOP, NOP) = MatrixXd::Identity(NOP, NOP);

	K.block(NOP, 2 * NOP, NOP, 1) = -1*MatrixXd::Ones(NOP,1);
	K.block(2 * NOP, NOP, 1, NOP) = R_c;

	K(2 * NOP, 2 * NOP) = L_f / dt;
	K(2 * NOP, 2 * NOP + 1) = -CalcR_dyn(t);
	K(2 * NOP + 1, 2 * NOP) = 1;
	K(2 * NOP + 1, 2 * NOP + 1) = 1;
}

void Circuit::CalcC(double t, int i_for_time) {		//Be careful. This fuction may use value out of index.
	double dt = TIME_TICK;

	C.setZero(NUMBER_OF_UNKOWNS);

	for (int i = 0; i < NOP; i++) {
		for (int j = 0; j < NOP; j++) {
			C(i) = M_ij(i, j) / dt * I_theta(j, i_for_time - 1);
		}
	}

	C(2 * NOP) = L_f / dt * I_L(i_for_time - 1);
	C(2 * NOP + 1) = CalcI_2(t);
}

//ダイナミック抵抗を計算（矩形波）
double Circuit::CalcR_dyn(double t) {
	double R_dyn;
	double i_2_period = 1 / i_2_freq;

	if (fmod(t, i_2_period) < (0.25*i_2_period - b_app_rate * 0.5*i_2_period))
		R_dyn = 0;
	else if (fmod(t, i_2_period) < (0.25*i_2_period + b_app_rate * 0.5*i_2_period))
		R_dyn = r_dyn_mag;
	else
		R_dyn = 0;

	return R_dyn;
}

//i2を計算（三角波）
double Circuit::CalcI_2(double t) {
	double i_2_period = 1 / i_2_freq;
	double pos = fmod(t, i_2_period);
	double a = i_2_mag / (0.25*i_2_period);

	if (pos < 0.25*i_2_period) {
		return a * pos;
	} else if (pos < 0.75 * i_2_period) {
		return -a * (pos - 0.5*i_2_period);
	} else {
		return a * (pos - i_2_period);
	}
}

void Circuit::CalcCircuit() {
	double dt = TIME_TICK;
	for (int i = 1; i < TIME_DIV_NUM; i++) {
		//t[i] = i * dt;

		//CalcC(i);
		//CalcM_ij();
		//CalcK(t[i]);

	}
}



