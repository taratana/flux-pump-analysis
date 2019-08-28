#pragma once
#include "headings.h"
#include "Eigen/Dense"
#include "Circuit.h"

using namespace std;
using namespace Eigen;


class DataFile {
private:
	ofstream output_file;
	int flag_for_style;

public:
	DataFile(string filename);
	~DataFile();
	void FirstOutput(string ss);

	void Output(double t, VectorXd i_theta, VectorXd i_r, double i_L, double i_B, double i_2, double r_dyn);
	void Output(MatrixXd &mat);
	void Output(VectorXd &vec);
	
};

