#include"DataFile.h"

using namespace std;
using namespace Eigen;

DataFile::DataFile(string filename) {
	output_file.open(filename);
}

DataFile::~DataFile() {
	output_file.close();
}

void DataFile::FirstOutput(string ss) {
	output_file << ss;
}

void DataFile::Output(double t, VectorXd  i_theta, VectorXd  i_r, double i_L, double i_B, double i_2, double r_dyn) {
	output_file << t << ",";
	for (int i = 0; i < i_theta.size(); i++)
		output_file << i_theta(i) << ",";
	for (int i = 0; i < i_r.size(); i++)
		output_file << i_r(i) << ",";
	output_file << i_L << ",";
	output_file << i_B << ",";
	output_file << i_2 << ",";
	output_file << r_dyn << endl;
}


//Output matrix data
void DataFile::Output(MatrixXd &mat){
	for (int i = 0; i < mat.rows(); i++) {
		for (int j = 0; j < mat.cols(); j++) {
			output_file << mat(i, j) << ",";
		}
		output_file << endl;
	}	
}

void DataFile::Output(VectorXd &vec) {
	for (int i = 0; i < vec.size(); i++) {
		output_file << vec(i) << endl;
	}
}


