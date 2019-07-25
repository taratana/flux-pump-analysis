#pragma once
#include<math.h>
#include<fstream>
#include<omp.h>
#include<iostream>
#include<iomanip>
#include<string>
#include<vector>
#include<Eigen/Dense>

static const double CRITICAL_CURRENT_OF_YBCO_TAPE = 95.0;
static const double NI_COIL_INDUCTANCE = 120e-6;
static const double FILTERING_COIL_INDUCTANCE = 16e-6;
static const double B_MAGNITUDE = 0.685;
static const double B_FRECUENCY = 28.0;
static const double i_2_MAGNITUDE = 120.0;
static const double i_2_FREQUENCY = 1/1.75; // = 0.5714
static const double R_c = 1.0e-6;
static const double R_JOINT = 176.5e-9;

static const double B_APP_LENGTH = 0.035;
static const double B_APP_WIDTH = 0.006;
static const double C_OF_DYNAMIC_RESISTENCE = 0.395;
static const double R_2 = 60e-6;	


static const double B_APP_RATE = 0.1;
static const double TIME_LENGTH = 3100;
static const double TIME_TICK = 0.001;
static const int TIME_DIV_NUM = (int)TIME_LENGTH / TIME_TICK;

static const int NOP = 12;	//Number of pancakes
static const int NUMBER_OF_UNKOWNS = 2 * NOP + 2;

static const int M_FOR_GMRES = 20;

static const double mu0 = 1.2566370614e-6;