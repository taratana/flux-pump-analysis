#include "circuit.h"

using namespace std;


circuit::circuit() {
	r_c = R_c;
	r_joint = R_JOINT;
	i_2_freq = i_2_FREQUENCY;
	i_2_mag = i_2_MAGNITUDE;
	l_coil = NI_COIL_INDUCTANCE;
	l_filter = FILTERING_COIL_INDUCTANCE;
	b_app_rate = B_APP_RATE;
	b_mag = B_MAGNITUDE;
	b_freq = B_FRECUENCY;
	b_app_length = B_APP_LENGTH;
	b_app_width = B_APP_WIDTH;
	r_dyn_c = C_OF_DYNAMIC_RESISTENCE;
	r_dyn_mag = 2 * b_app_length*b_app_width*b_freq*(b_mag + r_dyn_c * b_mag*b_mag) / CRITICAL_CURRENT_OF_YBCO_TAPE;
}

void circuit::setI2(double i_2_mag, double i_2_freq) {
	this->i_2_mag = i_2_mag;
	this->i_2_freq = i_2_freq;
}

void circuit::setp(double p) {
	this->b_app_rate = p;
}

void circuit::setRdyn(double rdyn) {
	this->r_dyn_mag = rdyn;
}

void circuit::CalcNoFilterCircuit(string output_filename) {
	double i_theta = 0, i_theta_tmp = 0;
	double i_r = 0;
	double k1, k2, k3, k4;
	std::ofstream output(output_filename);

	output << "#R_dyn:" << r_dyn_mag << endl;
	output << "#r_joint:" << r_joint << endl;
	output << "#r_c:" << r_c << endl;
	output << "#i_2_mag:" << i_2_mag << endl;
	output << "#i_2freq:" << i_2_freq << endl;
	output << "#l_coil:" << l_coil << endl;
	output << "#l_filter:" << l_filter << endl;
	output << "#b_app_rate:" << b_app_rate << endl;
	output << "#b_mag:" << b_mag << endl;
	output << "#b_app_length:" << b_app_length << endl;
	output << "#b_app_width:" << b_app_width << endl;
	output << "#r_dyn_c:" << r_dyn_c << endl;
	output << "#TIME_TICK:" << TIME_TICK << endl;
	output << "#Number of Time Step:" << (int)(TIME_LENGTH / TIME_TICK) << endl;
	output << "#t,i_theta,i_r,i_L,i_2,R_dyn" << endl;


	double t = 0;
	for (int i = 0; i <= (int)(TIME_LENGTH / TIME_TICK); i++) {
		k1 = TIME_TICK * FuncForNoFilterCircuit(t, i_theta);
		k2 = TIME_TICK * FuncForNoFilterCircuit(t + TIME_TICK * 0.5, i_theta + k1 * 0.5);
		k3 = TIME_TICK * FuncForNoFilterCircuit(t + TIME_TICK * 0.5, i_theta + k2 * 0.5);
		k4 = TIME_TICK * FuncForNoFilterCircuit(t + TIME_TICK, i_theta + k3);
		i_theta += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		t += TIME_TICK;

		i_r = l_coil / r_c * (i_theta - i_theta_tmp) / TIME_TICK;

		if (i % 4 == 0)
			output << t << "," << i_theta << "," << i_r << "," << i_theta + i_r << "," << CalcI_2(t) << "," << CalcR_dyn(t) << endl;
		cout << i << "Loop Finished." << endl;

		i_theta_tmp = i_theta;
	}
}

void circuit::CalcFilterExistCircuit(string output_filename) {
	double i_theta = 0, i_theta_d = 0, i_theta_tmp = 0;
	double i_r;
	double k1, k2, k3, k4;
	double l1, l2, l3, l4;
	std::ofstream output(output_filename);

	output << "#R_dyn:" << r_dyn_mag << endl;
	output << "#r_joint:" << r_joint << endl;
	output << "#r_c:" << r_c << endl;
	output << "#i_2_mag:" << i_2_mag << endl;
	output << "#i_2_freq:" << i_2_freq << endl;
	output << "#l_coil:" << l_coil << endl;
	output << "#l_filter:" << l_filter << endl;
	output << "#b_app_rate:" << b_app_rate << endl;
	output << "#b_mag:" << b_mag << endl;
	output << "#b_app_length:" << b_app_length << endl;
	output << "#b_app_width:" << b_app_width << endl;
	output << "#r_dyn_c:" << r_dyn_c << endl;
	output << "#TIME_TICK:" << TIME_TICK << endl;
	output << "#Number of Time Step:" << (int)(TIME_LENGTH / TIME_TICK) << endl;
	output << "#t,i_theta,i_r,i_L,i_2,R_dyn" << endl;

	double t = 0;
	for (int i = 0; i <= (int)(TIME_LENGTH / TIME_TICK); i++) {
		k1 = TIME_TICK * FuncForFilterCircuit1(t, i_theta, i_theta_d);
		l1 = TIME_TICK * FuncForFilterCircuit2(t, i_theta, i_theta_d);
		k2 = TIME_TICK * FuncForFilterCircuit1(t + TIME_TICK * 0.5, i_theta + k1 * 0.5, i_theta_d + l1 * 0.5);
		l2 = TIME_TICK * FuncForFilterCircuit2(t + TIME_TICK * 0.5, i_theta + k1 * 0.5, i_theta_d + l1 * 0.5);
		k3 = TIME_TICK * FuncForFilterCircuit1(t + TIME_TICK * 0.5, i_theta + k2 * 0.5, i_theta_d + l2 * 0.5);
		l3 = TIME_TICK * FuncForFilterCircuit2(t + TIME_TICK * 0.5, i_theta + k2 * 0.5, i_theta_d + l2 * 0.5);
		k4 = TIME_TICK * FuncForFilterCircuit1(t + TIME_TICK, i_theta + k3, i_theta_d + l3);
		l4 = TIME_TICK * FuncForFilterCircuit2(t + TIME_TICK, i_theta + k3, i_theta_d + l3);
		i_theta += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		i_theta_d += (l1 + 2 * l2 + 2 * l3 + l4) / 6;
		t += TIME_TICK;

		i_r = l_coil / r_c * (i_theta - i_theta_tmp) / TIME_TICK;

		if (i % 4 == 0)
			output << t << "," << i_theta << "," << i_r << "," << i_theta + i_r << "," << CalcI_2(t) << "," << CalcR_dyn(t) << endl;
		cout << i << "Loop Finished." << endl;

		i_theta_tmp = i_theta;


	}
}

//ダイナミック抵抗を計算（矩形波）
double circuit::CalcR_dyn(double t) {
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
double circuit::CalcI_2(double t) {
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

double circuit::FuncForNoFilterCircuit(double t, double i_theta) {
	double a, b;
	a = (r_c + r_joint + CalcR_dyn(t)) / r_c * l_coil;
	b = r_joint + CalcR_dyn(t);
	return CalcR_dyn(t) / a * CalcI_2(t) - b / a * i_theta;
}

double circuit::FuncForFilterCircuit1(double t, double i_theta, double i_theta_d) {
	return i_theta_d;
}

double circuit::FuncForFilterCircuit2(double t, double i_theta, double i_theta_d) {
	double a, b, c;
	a = (r_c + r_joint + CalcR_dyn(t)) / r_c * l_coil;
	b = r_joint + CalcR_dyn(t);
	c = l_filter * l_coil / r_c;

	return CalcR_dyn(t) / c * CalcI_2(t) - (a + l_filter) / c * i_theta_d - b / c * i_theta;
}

