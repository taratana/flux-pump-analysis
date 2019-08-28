#include "Circuit.h"

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

	nop = NOP;
	nou = NUMBER_OF_UNKNOWNS;
}

//Initializing arrays
void Circuit::InitMatrices() {
	cout << "Initializing matrices." << endl;
	time = VectorXd::Zero(TIME_DIV_NUM);
	I_L = VectorXd::Zero(TIME_DIV_NUM);
	I_B = VectorXd::Zero(TIME_DIV_NUM);
	R_c = VectorXd::Zero(nop);
	I_theta = MatrixXd::Zero(nop, TIME_DIV_NUM);
	I_r = MatrixXd::Zero(nop, TIME_DIV_NUM);

	C = VectorXd::Zero(nou);
	K = MatrixXd::Zero(nou, nou);
	M_ij = MatrixXd::Zero(nop, nop);
}

void Circuit::setI2(double i_2_mag, double i_2_freq) {
	this->i_2_mag = i_2_mag;
	this->i_2_freq = i_2_freq;
}

void Circuit::setp(double p) {
	this->b_app_rate = p;
}

void Circuit::setFilterInductance(double L_f) {
	this->L_f = L_f;
}

void Circuit::setR_dyn(double r_dyn) {
	this->r_dyn_mag = r_dyn;
}

void Circuit::setR_c(double r_c) {
	for (int i = 0; i < nop; i++)R_c(i) = r_c;
}

void Circuit::setResistance(double R_c, double R_joint) {

}

void Circuit::setNOP(int nop) {
	this->nop = nop;
	this->nou = 2 * nop + 2;
}

void Circuit::getConvergenceInfo(double * conv_t, double * conv_v) {
	*conv_t = this->convergence_time;
	*conv_v = this->convergence_value;
}

void Circuit::CalcM_ij() {
	string line;
	ifstream in;
#ifdef LOAD_INDUCTANCE_MATRIX
	in.open("ExpData/Inductance/Inductance-Matrix" + to_string(nop) + ".csv");
#endif // LOAD_INDUCTANCE_MATRIX

	//If csv file exist, the file is loaded.
	if (in.is_open()) {
		for (int i = 0; i < nop; i++) {
			getline(in, line);
			string token;
			istringstream stream(line);

			for (int j = 0; j < nop; j++) {
				getline(stream, token, ',');
				M_ij(i, j) = atof(token.c_str());
			}
		}
	} else {
		int i, j;
		double r_w = REBCO_WIDTH;
		double i_w = INSULATOR_WIDTH;

		cout << "Inductance matrix calculation started." << endl;
		for (j = 0; j < nop; j++) {
			inductance.sc1->setCoilParameter(INNNER_DIAMETER / 2, OUTER_DIAMETER / 2, r_w*j + i_w * j, r_w*(j + 1) + i_w * j, NUMBER_OF_TURNS);
			for (i = 0; i < nop; i++) {
				inductance.sc2->setCoilParameter(INNNER_DIAMETER / 2, OUTER_DIAMETER / 2, r_w*i + i_w * i, r_w*(i + 1) + i_w * i, NUMBER_OF_TURNS);
				M_ij(i, j) = inductance.M();
			}
			std::cout << "Coil" << j + 1 << " done." << std::endl;
		}
#ifndef DONT_OUTPUT_CURRENTS_DATA
		DataFile *opd = new DataFile("ExpData/Inductance/Inductance-Matrix" + to_string(nop) + ".csv");
		opd->Output(M_ij);
		delete opd;
#endif // DONT_OUTPUT_CURRENTS_DATA

	}
	total_inductance = M_ij.sum();
}

void Circuit::CalcK(double t) {
	double dt = TIME_TICK;

	K.setZero(nou, nou);

	//cout << "Calculating Coefficient matrix." << endl;

	K.block(0, 0, nop, nop) = M_ij / dt;
	K.block(nop, 0, nop, nop) = MatrixXd::Identity(nop, nop);
	K.block(0, nop, nop, nop) = -1 * R_c.asDiagonal();
	K.block(nop, nop, nop, nop) = MatrixXd::Identity(nop, nop);

	for (int i = 0; i < nop; i++) {
		K(2 * nop, i) = R_c(i);
	}
	K.block(nop, 2 * nop, nop, 1) = -1 * MatrixXd::Ones(nop, 1);
	K(2 * nop, 2 * nop) = -(L_f / dt + R_c.sum());
	K(2 * nop, 2 * nop + 1) = CalcR_dyn(t);
	K(2 * nop + 1, 2 * nop) = 1;
	K(2 * nop + 1, 2 * nop + 1) = 1;
	

#ifndef DONT_OUTPUT_MISCELLANEOUS_DATA
	DataFile *opd = new DataFile("ExpData/miscellaneous/Coefficient-Matrix.csv");
	opd->Output(K);
	delete opd;
#endif // DONT_OUTPUT_MISCELLANEOUS_DATA
}

void Circuit::CalcC(double t, int i_at_t) {
	//Be careful. This fuction may use value out of index at i_for_time=0.
	double dt = TIME_TICK;

	C.setZero();

	//cout << "Calculating C." << endl;

	for (int i = 0; i < nop; i++) {
		for (int j = 0; j < nop; j++) {
			C(i) += M_ij(i, j) / dt * I_theta(j, i_at_t - 1);
		}
	}

	C(2 * nop) = -L_f / dt * I_L(i_at_t - 1);
	C(2 * nop + 1) = CalcI_2(t);


#ifndef DONT_OUTPUT_MISCELLANEOUS_DATA
	DataFile *opd = new DataFile("ExpData/miscellaneous/C.csv");
	opd->Output(C);
	delete opd;
#endif // DONT_OUTPUT_MISCELLANEOUS_DATA

}

void Circuit::CalcRc() {
	double Rc_of_one_coil = 0;

	cout << "Calculating Rc." << endl;
	for (int i = 1; i <= NUMBER_OF_TURNS - 1; i++) {
		Rc_of_one_coil += 1 / (INNNER_DIAMETER + (i + 0.5)*REBCO_THICKNESS);
	}

	Rc_of_one_coil *= CONTACT_RESISTIVITY / (2 * M_PI*REBCO_WIDTH);

	for (int i = 0; i < R_c.size(); i++) {
		R_c(i) = Rc_of_one_coil;
	}

#ifndef DONT_OUTPUT_MISCELLANEOUS_DATA
	DataFile *opd = new DataFile("ExpData/miscellaneous/Contact-Resistive.csv");
	opd->Output(R_c);
	delete opd;
#endif // DONT_OUTPUT_MISCELLANEOUS_DATA
}

//ダイナミック抵抗を計算（矩形波）
double Circuit::CalcR_dyn(double t) {
	double R_dyn;
	double i_2_period = 1 / Circuit::i_2_freq;

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


void Circuit::CalcCircuit(string filename) {
#ifndef DONT_OUTPUT_CURRENTS_DATA
	opd = new DataFile(filename + ".csv");
	opd->FirstOutput("#t,");
	for (int i = 0; i < nop; i++)
		opd->FirstOutput("I_theta" + to_string(i) + ",");
	for (int i = 0; i < nop; i++)
		opd->FirstOutput("I_r" + to_string(i) + ",");
	opd->FirstOutput("I_L,I_B,I_2,R_dyn\n");
#endif // DONT_OUTPUT_CURRENTS_DATA
	
	InitMatrices();
	CalcRc();
	CalcM_ij();

	double dt = TIME_TICK;
	VectorXd I = VectorXd::Ones(nou);
	for (int i = 1; i < TIME_DIV_NUM; i++) {
		time(i) = (double)i * dt;
		cout << "\n<LOOP i=" << i << " started (t=" << time(i) << ")>" << endl;
		CalcC(time(i), i);
		CalcK(time(i));
		//cout << "K" << endl << K << endl;
		//cout << "C" << endl << C << endl;
		//I.setRandom();
		//I.setLinSpaced(0.1, 10);
		//I.setOnes();

		I = GMRES(K, C, I, EPS);

		for (int j = 0; j < nop; j++) {
			I_theta(j, i) = I(j);
			I_r(j, i) = I(j + nop);
		}
		I_L(i) = I(2 * nop);
		I_B(i) = I(2 * nop + 1);

#ifndef DONT_OUTPUT_CURRENTS_DATA
		if (i % 10 == 0)
			opd->Output(time(i), I_theta.col(i), I_r.col(i), I_L(i), I_B(i), CalcI_2(time(i)), CalcR_dyn(time(i)));
#endif // DONT_OUTPUT_CURRENTS_DATA

	}

	//Convergence Calculation
	CalcConvergenceInfo();

#ifndef DONT_OUTPUT_CURRENTS_DATA
	PrintInfo(filename);
	delete opd;
#endif // DONT_OUTPUT_CURRENTS_DATA

}

void Circuit::CalcConvergenceInfo() {
	int maximum_pancake = 0;
	convergence_value = I_theta(0, TIME_DIV_NUM - 1);
	for (int i = 0; i < nop; i++) {
		if (I_theta(i, TIME_DIV_NUM - 1) > convergence_value) {
			convergence_value = I_theta(i, TIME_DIV_NUM - 1);
			maximum_pancake = i;
		}
	}
	for (int i = 0; i < TIME_DIV_NUM; i++) {
		if (I_theta(maximum_pancake, i) > 0.95*convergence_value) {
			convergence_time = i * TIME_TICK;
			break;
		}
	}
}

void Circuit::PrintInfo(string filename) {
	ofstream info_output(filename + ".txt");

	info_output << "#[Simulation Results]" << endl;
	info_output << "#CONVERGENCE_VALUE:" << convergence_value << endl;
	info_output << "#CONVERGENCE_TIME:" << convergence_time << endl;
	info_output << endl;

	info_output << "#[Simulation Condituions, etc]" << endl;
	info_output << "#NUMBER_OF_UNKNOWNS:" << to_string(nou) << endl;
	info_output << "#TIME_LENGTH:" << TIME_LENGTH << endl;
	info_output << "#TIME_TICK:" << TIME_TICK << endl;
	info_output << "#TIME_DIV_NUM:" << TIME_DIV_NUM << endl;
	info_output << "#GMRES_EPS:" << EPS << endl;
	info_output << "#GMRES_M:" << M_FOR_GMRES << endl;
	info_output << endl;

	info_output << "#[Pancake Coil Parameters]" << endl;
	info_output << "#NOP:" << nop << endl;
	info_output << "#REBCO_WIDTH:" << REBCO_WIDTH << endl;
	info_output << "#REBCO_THICKNESS:" << REBCO_THICKNESS << endl;
	info_output << "#REBCO_INNNER_DIAMETER:" << INNNER_DIAMETER << endl;
	info_output << "#REBCO_OUTER_DIAMETER:" << OUTER_DIAMETER << endl;
	info_output << "#NUMBER_OF_TURNS:" << NUMBER_OF_TURNS << endl;
	info_output << "#CONTACT_RESISTIVE:" << R_c(0) << endl;
	info_output << "#CONTACT_RESISTIVITY:" << CONTACT_RESISTIVITY << endl;
	info_output << "#INSULATOR_WIDTH:" << INSULATOR_WIDTH << endl;
	info_output << endl;

	info_output << "#[Flux Pump Parameters]" << endl;
	info_output << "#FILTER_INDUCTANCE:" << L_f << endl;
	info_output << "#TOTAL_INDUCTANCE_OF_PANCAKES:" << total_inductance << endl;
	info_output << "#R_dyn:" << r_dyn_mag << endl;
	info_output << "#I2_MAGNITUDE:" << i_2_mag << endl;
	info_output << "#I2_FREQUENCY:" << i_2_freq << endl;
	info_output << "#B_MAGNITUDE:" << b_mag << endl;
	info_output << "#B_FREQUENCY:" << b_freq << endl;
	info_output << "#B_APP_RATE:" << B_APP_RATE << endl;
	info_output << "#B_APP_LENGTH:" << b_app_length << endl;
	info_output << "#B_APP_WIDTH:" << b_app_width << endl;
	info_output << "#C_OF_DYNAMIC_RESISTANCE:" << r_dyn_c << endl;
	info_output << "#REBCO_TAPE_CRITICAL_CURRENT:" << CRITICAL_CURRENT_OF_YBCO_TAPE << endl;

	info_output.close();
}

