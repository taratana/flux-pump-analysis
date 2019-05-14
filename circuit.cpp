#include "circuit.h"

using namespace std;


void circuit::CalcNoFilterCircuit() {
	double i_theta = 0;
	double k1, k2, k3, k4;
	std::ofstream output("output-without-filter.csv");

	double t = 0;
	for(int i=0;i<=(int)(TIME_LENGTH/TIME_TICK);i++){
		k1 = TIME_TICK * FuncForNoFilterCircuit(t, i_theta);
		k2 = TIME_TICK * FuncForNoFilterCircuit(t + TIME_TICK * 0.5, i_theta + k1 * 0.5);
		k3 = TIME_TICK * FuncForNoFilterCircuit(t + TIME_TICK * 0.5, i_theta + k2 * 0.5);
		k4 = TIME_TICK * FuncForNoFilterCircuit(t + TIME_TICK, i_theta + k3);
		i_theta += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		t += TIME_TICK;

		if (i % 3 == 0)
			output << t << "," << i_theta << "," << CalcI_2(t) << "," << CalcR_dyn(t) << endl;
		cout  << i << "Loop Finished." << endl;
		
	}
}

void circuit::CalcFilterExistCircuit() {
	double i_theta = 0, i_theta_d = 0;
	double k1, k2, k3, k4;
	double l1, l2, l3, l4;
	std::ofstream output("output-with-filter.csv");



	double t = 0;
	for (int i = 0; i <= (int)(TIME_LENGTH / TIME_TICK); i++) {
		k1 = TIME_TICK * FuncForFilterCircuit1(t, i_theta,i_theta_d);
		l1 = TIME_TICK * FuncForFilterCircuit2(t, i_theta,i_theta_d);
		k2 = TIME_TICK * FuncForFilterCircuit1(t + TIME_TICK * 0.5, i_theta + k1 * 0.5, i_theta_d + l1 * 0.5);
		l2 = TIME_TICK * FuncForFilterCircuit2(t + TIME_TICK * 0.5, i_theta + k1 * 0.5, i_theta_d + l1 * 0.5);
		k3 = TIME_TICK * FuncForFilterCircuit1(t + TIME_TICK * 0.5, i_theta + k2 * 0.5, i_theta_d + l2 * 0.5);
		l3 = TIME_TICK * FuncForFilterCircuit2(t + TIME_TICK * 0.5, i_theta + k2 * 0.5, i_theta_d + l2 * 0.5);
		k4 = TIME_TICK * FuncForFilterCircuit1(t + TIME_TICK, i_theta + k3, i_theta_d + l3);
		l4 = TIME_TICK * FuncForFilterCircuit2(t + TIME_TICK, i_theta + k3, i_theta_d + l3);


		i_theta += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		i_theta_d += (l1 + 2 * l2 + 2 * l3 + l4) / 6;
		t += TIME_TICK;

		if(i%3==0)
			output << t << "," << i_theta << "," << i_theta_d<<"," << CalcI_2(t) << "," << CalcR_dyn(t)  << endl;
		cout << i << "Loop Finished." << endl;

	}
}

//ダイナミック抵抗を計算（矩形波）
double circuit::CalcR_dyn(double t) {
	double R_dyn;
	double i_2_period = 1 / i_2_FREQUENCY;

	if(fmod(t,i_2_period)<(0.25*i_2_period-B_APP_RATE*0.5*i_2_period))
		R_dyn = 0;
	else if(fmod(t, i_2_period) < (0.25*i_2_period + B_APP_RATE*0.5*i_2_period))
		R_dyn = 2 * B_APP_LENGTH*B_APP_WIDTH*B_FRECUENCY*B_MAGNITUDE*B_MAGNITUDE / CRITICAL_CURRENT_OF_YBCO_TAPE;
	else
		R_dyn = 0;

	return R_dyn;
}

//i2を計算（三角波）
double circuit::CalcI_2(double t) {
	double i_2_period = 1 / i_2_FREQUENCY;
	double pos = fmod(t,i_2_period);
	double a = i_2_MAGNITUDE / (0.25*i_2_period);

	if (pos < 0.25*i_2_period) {
		return a * pos;
	} else if (pos < 0.75 * i_2_period) {
		return -a * (pos - 0.5*i_2_period);
	} else {
		return a * (pos - i_2_period);
	}
}

double circuit::FuncForNoFilterCircuit(double t, double i_theta) {
	double a, b;
	a = (R_c + R_JOINT + CalcR_dyn(t)) / R_c * NI_COIL_INDUCTANCE;
	b = R_JOINT + CalcR_dyn(t);
	return CalcR_dyn(t) / a * CalcI_2(t) - b / a * i_theta;
}

double circuit::FuncForFilterCircuit1(double t, double i_theta,double i_theta_d) {
	return i_theta_d;
}

double circuit::FuncForFilterCircuit2(double t, double i_theta, double i_theta_d) {
	double a, b, c;
	a = (R_c + R_JOINT + CalcR_dyn(t)) / R_c*NI_COIL_INDUCTANCE;
	b = R_JOINT + CalcR_dyn(t);
	c = FILTERING_COIL_INDUCTANCE * NI_COIL_INDUCTANCE / R_c;

	return CalcR_dyn(t) / c * CalcI_2(t) - (a + FILTERING_COIL_INDUCTANCE) / c * i_theta_d - b / c * i_theta;
}

