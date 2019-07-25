#include"OutputData.h"

using namespace std;



OutputData::OutputData(string filename, EOutputStyle output_style) {
	flag_for_style = output_style;
	output_file.open(filename);
	FirstOutput();
}

OutputData::~OutputData() {
	output_file.close();
}

void OutputData::FirstOutput() {
	output_file << "#NOP=" << NOP << endl;
	//‚±‚±‚É’Ç‰Á‚·‚é•K—v‚ ‚èD

	switch (flag_for_style) {
	case ALL_CURRENT:
		output_file << "#For Gnuplot" << endl;
		output_file << "#t,I_theta1,...,IthetaNOP,Ir1,...,IrNOP,IB,IL" << endl;
		break;
	default:
		break;
	}
}

void OutputData::outputData(double a,double b,double c) {
	output_file << a << "," << b << "," << c << endl;
}
