#include "Inductance.h"


Inductance::SingleCoil::SingleCoil() {
	setDivNum(N0, NR, NZ);
}

double Inductance::SingleCoil::Hr(){
	if (Rp == 0.0) {
		return 0.0;
	}else {
		/*return deployFunction(a1, a2, z1, z2, Hr1)+deployFunction(a1,a2,z1,z2,Hr2);*/
		return Hr1(a2, z2) - Hr1(a1, z2) - Hr1(a2, z1) + Hr1(a1, z1)
			+ Hr2(a2, z2) - Hr2(a1, z2) - Hr2(a2, z1) + Hr2(a1, z1);
	}
}

double Inductance::SingleCoil::Hz(){
	if (Rp == 0.0) {
		/*return deployFunction(a1, a2, z1, z2, Hz0);*/
		return Hz0(a2, z2) - Hz0(a1, z2) - Hz0(a2, z1) + Hz0(a1, z1);
	}else {
		/*return deployFunction(a1, a2, z1, z2, Hz1) + deployFunction(a1, a2, z1, z2, Hz2) + deployFunction(a1, a2, z1, z2, Hz3);*/
		return Hz1(a2, z2) - Hz1(a1, z2) - Hz1(a2, z1) + Hz1(a1, z1)
			+ Hz2(a2, z2) - Hz2(a1, z2) - Hz2(a2, z1) + Hz2(a1, z1)
			+ Hz3(a2, z2) - Hz3(a1, z2) - Hz3(a2, z1) + Hz3(a1, z1);
	}
}

double Inductance::SingleCoil::Bz(double Rp, double Zp) {
	setRpZp(Rp, Zp);
	return mu0 * Hz();
}

double Inductance::SingleCoil::Hz0(double a, double z) {
	return J / 2 * (z - Zp) * (2 * log(10) + log(sqrt(a * a + (z - Zp) * (z - Zp)) + a));
}

/*c++ではローカル関数は禁止されているので、
ラムダ式を使う*/
double Inductance::SingleCoil::Hr1(double a, double z) {
	auto integrand = [this, a, z](double theta)->double {
		return R(a, z, theta) * cos(theta);
	};
	return 1 / (4 * M_PI) * J * gaussIntegration.integrate(0.0, 2 * M_PI, integrand);
}

double Inductance::SingleCoil::Hr2(double a, double z){
	if (Zp != z) {
		return Hr2_3(a, z);
	}
	else if (Rp < a) {
		return Hr2_3(a, z);
	}
	else if (Rp > a) {
		return Hr2_2(a, z);
	}
	else if(Rp == a){
		return Hr2_1(a, z);
	}
}

//Zp=z, Rp=a
double Inductance::SingleCoil::Hr2_1(double a, double z) {
	auto integrand = [this, a, z](double theta)->double {
		return cos(theta) * cos(theta) * log(sqrt(2) + sqrt(1 - cos(theta)));
	};
	return J / (4 * M_PI) * Rp * (
		M_PI * log(Rp) - M_PI / 4 * (1 + 2 * log(2)) + 2 * M_PI * log(10) +
		gaussIntegration.integrate(0.0, 2 * M_PI, integrand)
		);
}

//Zp=z, Rp>a
double Inductance::SingleCoil::Hr2_2(double a, double z) {
	auto integrand1 = [this, a, z](double theta)->double {
		return cos(theta) * cos(theta) *
			log(R(a, z, theta) - a + Rp * cos(theta));
	};
	auto integrand2 = [this, a, z](double theta)->double {
		return cos(theta) * cos(theta) *
			log(R(a, z, theta) + a - Rp * cos(theta));
	};
	return J / (2 * M_PI) * Rp * (M_PI / 2 * log(Rp) - M_PI / 4 * (1 + 2 * log(2))
		+ M_PI * log(10) - gaussIntegration.integrate(0.0, M_PI / 2, integrand1) + gaussIntegration.integrate(M_PI / 2, M_PI, integrand2));
}

//Zp=z, Rp<a or Zp!=z
double Inductance::SingleCoil::Hr2_3(double a, double z) {
	auto integrand = [this, a, z](double theta)->double {
		return cos(theta) * cos(theta) * log(
			R(a, z, theta) + a - Rp * cos(theta)
		);
	};
	return J / (4 * M_PI) * Rp * (2 * M_PI * log(10) + gaussIntegration.integrate(0.0, M_PI, integrand));
}

double Inductance::SingleCoil::Hz1(double a, double z) {
	if (Zp != z) {
		return Hz1_2(a, z);
	} else if(Zp == z){
		return Hz1_1(a, z);
	}
}

//Zp=0
double Inductance::SingleCoil::Hz1_1(double a, double z) {
	return 0.0;
}

//Zp!=0
double Inductance::SingleCoil::Hz1_2(double a, double z) {
	auto integrand = [this, a, z](double theta)->double {
		return log(R(a, z, theta) + a - Rp * cos(theta));
	};
	return J / (4 * M_PI) * (z - Zp) * (4 * M_PI * log(10) + gaussIntegration.integrate(0.0, 2 * M_PI, integrand));
}

double Inductance::SingleCoil::Hz2(double a, double z) {
	if (Zp == z) {
		return Hz2_1(a, z);
	} else if (Rp == a) {
		return Hz2_2(a, z);
	} else if(Rp != a){
		return Hz2_3(a, z);
	}
}

//Zp=z
double Inductance::SingleCoil::Hz2_1(double a, double z) {
	return 0.0;
}

//Zp!=z, Rp=a
double Inductance::SingleCoil::Hz2_2(double a, double z) {
	auto integrand = [this, a, z](double theta)->double {
		return cos(theta) * log(R(a, z, theta) + fabs(z - Zp));
	};
	return -J / (4 * M_PI) * Rp * (z - Zp) / fabs(z - Zp) * (M_PI + gaussIntegration.integrate(0.0, 2 * M_PI, integrand));

}

//Zp!=z, Rp!=a
double Inductance::SingleCoil::Hz2_3(double a, double z) {
	auto integrand = [this, a, z](double theta)->double {
		return cos(theta) * log((R(a, z, theta) - fabs(z - Zp)) /(R(a, z, theta) + fabs(z - Zp)));
	};
	return J / (8 * M_PI) * Rp * (z - Zp) / fabs(z - Zp) * gaussIntegration.integrate(0.0, 2 * M_PI, integrand);
}


double Inductance::SingleCoil::Hz3(double a, double z) {
	auto g = [this, a, z](double theta)->double {
		if (theta == 0 || theta == M_PI || theta == 2 * M_PI) {
			return 0.0;
		}
		else {
			return sin(theta) * atan((a - Rp * cos(theta)) * (z - Zp) / (Rp * R(a, z, theta) * sin(theta)));
		}
	};
	return -J / (4 * M_PI) * Rp * gaussIntegration.integrate(0.0, 2 * M_PI, g);
}


double Inductance::SingleCoil::R(double a, double z, double theta) {
	return sqrt(a * a - 2 * a * Rp * cos(theta) + Rp * Rp + (z - Zp) * (z - Zp));
}


double Inductance::SingleCoil::CalcJ(){
	return number_of_turn * 1.0 / ((a2 - a1) * (z2 - z1));
}

void Inductance::SingleCoil::setRpZp(double _Rp,double _Zp){
	Rp = _Rp;
	Zp = _Zp;
}

void Inductance::SingleCoil::setDivNum(int n0, int nr, int nz) {
	this->n0 = n0;
	this->nr = nr;
	this->nz = nz;
}

void Inductance::SingleCoil::setCoilParameter(double _a1, double _a2, double _z1, double _z2, int turn) {
	a1 = _a1;
	a2 = _a2;
	z1 = _z1;
	z2 = _z2;

	number_of_turn = turn;

	J = CalcJ();
}


Inductance::Inductance() {
	sc1 = new SingleCoil();
	sc2 = new SingleCoil();
}


void Inductance::CalcInductance(Eigen::MatrixXd & M) {
	//int i, j;
	//double r_w = REBCO_WIDTH;
	//double i_w = INSULATOR_WIDTH;
	//for (j = 0; j < M.rows(); j++) {
	//	sc1->setCoilParameter(INNNER_DIAMETER / 2, OUTER_DIAMETER / 2, r_w*j + i_w * j, r_w*(j + 1) + i_w * j, NUMBER_OF_TURN);
	//	for (i = 0; i < M.cols(); i++) {
	//		sc2->setCoilParameter(INNNER_DIAMETER / 2, OUTER_DIAMETER / 2, r_w*i + i_w * i, r_w*(i + 1) + i_w * i, NUMBER_OF_TURN);
	//		M(i, j) = this->M();
	//	}
	//	std::cout << "coil" << j+1 << " done." << std::endl;
	//}
}





double Inductance::M() {
	int h, i, j;
	double temp = 0, temp1, temp2;

	for (h = 1; h <= sc2->nz; h++) {
		temp1 = 0.0;
		for (i = 1; i <= sc2->n0; i++) {
			temp1 += ri(i) * sc1->Bz(ri(i), zh(h));;
		}
		temp1 *= sc2-> a1 / sc2->n0;


		temp2 = 0.0;
		for (j = 1; j <= sc2->nr; j++) {
			temp2 += rj(j) * (sc2->a2 - rj(j)) * sc1->Bz(rj(j), zh(h));
		}
		temp2 *= 1.0 / sc2->nr;

		temp += temp1 + temp2;
		//std::cout << h << std::endl;
	}
	temp *= 2 * M_PI * sc2->number_of_turn / sc2->nz;

	return temp;
}

double Inductance::ri(int i) {
	return (2.0 * i - 1) * sc2->a1 / (2.0 * sc2->n0);
}

double Inductance::rj(int j) {
	return sc2->a1 + (2.0 * j - 1) * (sc2->a2 - sc2->a1) / (2.0 * sc2->nr);
}

double Inductance::zh(int h) {
	return sc2->z1 + (2.0 * h - 1) * (sc2->z2 - sc2->z1) / (2.0 * sc2->nz);
}