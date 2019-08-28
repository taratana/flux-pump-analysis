#pragma once
#define _USE_MATH_DEFINES
//#define EIGEN_DONT_PARALLELIZE
#include<math.h>
#include<omp.h>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<string>
#include "Eigen/Dense"

#define LOAD_INDUCTANCE_MATRIX
//#define DONT_OUTPUT_CURRENTS_DATA
//#define DONT_OUTPUT_MISCELLANEOUS_DATA

//Flux pump parameters
static const double CRITICAL_CURRENT_OF_YBCO_TAPE = 95.0;
static const double FILTERING_COIL_INDUCTANCE = 160e-6;
static const double B_MAGNITUDE = 0.685;
static const double B_FRECUENCY = 25.0;//28.0
static const double i_2_MAGNITUDE = 120.0;
static const double i_2_FREQUENCY = 1/2.0; // = 0.5

static const double B_APP_LENGTH = 0.100;	//0.065
static const double B_APP_WIDTH = 0.010;	//0.001
static const double C_OF_DYNAMIC_RESISTENCE = 0.395;

static const double B_APP_RATE = 0.1;

//Pancake parameters
static const int NOP = 3;	//Number of pancakes
static const int NUMBER_OF_TURNS = 100;
static const int NUMBER_OF_UNKNOWNS = 2 * NOP + 2;
static const double REBCO_WIDTH = 4e-3;
static const double INSULATOR_WIDTH = 1e-3;
static const double INNNER_DIAMETER = 60e-3;
static const double OUTER_DIAMETER = 79.2e-3;
static const double REBCO_THICKNESS = 0.096e-3;
static const double CONTACT_RESISTIVITY = 70e-10; // [Ohm m^2] -> 70[u Ohm cm^2]

//Constants
static const double mu0 = 1.2566370614e-6;

//Simulation conditions
static const double TIME_LENGTH = 30000;
static const double TIME_TICK = 0.01;
static const int TIME_DIV_NUM = (int)TIME_LENGTH / TIME_TICK;
static const double EPS = 1e-8;
static const int DATA_REDUCTION_RATE = 15;