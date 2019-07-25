#include"headings.h"

using namespace std;

enum EOutputStyle {
	ALL_CURRENT,

};

class OutputData {
private:
	ofstream output_file;
	int flag_for_style;

public:
	OutputData(string filename,EOutputStyle output_style);
	~OutputData();
	void FirstOutput();
	void outputData(double a,double b,double c);

};

