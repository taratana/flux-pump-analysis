#include<iostream>
#include"circuit.h"

using namespace std;

class StopWatch {
private:
	time_t  st_t, fin_t;  /*���Ԍv���p�i�����ԁj */
	clock_t st_ct, fin_ct; /*���Ԍv���p�iCPU���ԁj*/

public:
	void start() { /*���Ԍv���J�n*/
		st_ct = clock();
		time(&st_t);
	}
	void stop() { /*���Ԍv���I��*/
		fin_ct = clock();
		time(&fin_t);
	}
	double getTime() {
		return difftime(fin_t, st_t);
	}
};


int main() {
	StopWatch record_time;
	record_time.start();

	circuit c;
	double r_dyn;
	for(int i =0;i<5;i++){
		if (i == 0)
			r_dyn = 10e-6;
		else if (i == 1)
			r_dyn = 50e-6;
		else if (i == 2)
			r_dyn = 100e-6;
		else if (i == 3)
			r_dyn = 500e-6;
		else if (i == 4)
			r_dyn = 1000e-6;

		char output_name[20];
		sprintf(output_name, "with-filter-rdyn-%d.csv",i+1);
		c.setRdyn(r_dyn);
		c.CalcFilterExistCircuit(output_name);
	}


	record_time.stop();
	double elapsed = record_time.getTime();
	cout << "Time Elapsed:" << endl << (int)elapsed / 3600 << "hrs " << ((int)elapsed % 3600) / 60.0 << "min" << endl;
	
}