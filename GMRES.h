/*
This program solve simultaneous equation by GMRES(m).
But in iterations, EIGEN "normalized()" function may return same vector
because the norm is small enough. So normalization is done primitively.
You can use this solver like belelow:
	x = GMRES(A,b,x0,1e-8)
This program is made by Mato in 2019.
*/


#pragma once
#include<iostream>
#include<iomanip>

using namespace std;



static const int M_FOR_GMRES = NUMBER_OF_UNKNOWNS;
static const int KMAX_FOR_GMRES = 3000;
static const double ZERO_CHECK = 1e-30;
//#define DEBUG_GMRES1
//#define DEBUG_GMRES2
#define DIAGONAL_NORMALIZATION




// GMRES
template<class Matrix, class Vector>
Vector GMRES(Matrix &A, Vector &b, Vector &x0, double tol, int kmax, int m) {
	if ((A.cols() != A.rows()) || (b.size() != x0.size()) || (A.cols() != x0.size())) {
		std::cout << "Matrix, vector size error." << endl;
		exit(EXIT_FAILURE);	
	}

	//Get dimension of matrices.
	int dim = x0.size();

	double temp_norm = 0;
	int loop_count = 0;
	int Av_same_as_v = 0;	//Flag for cheacking Av=v
	Matrix h = MatrixXd::Zero(m + 1, m);
	Matrix r = MatrixXd::Zero(m, m);
	Matrix v = MatrixXd::Zero(dim, m + 1);
	Vector r0 = VectorXd::Zero(dim);
	Vector Av = VectorXd::Zero(dim);
	Vector temp = VectorXd::Zero(dim);
	Vector xm = VectorXd::Zero(dim);
	Vector g = VectorXd::Zero(m + 1);
	Vector c = VectorXd::Zero(m);
	Vector s = VectorXd::Zero(m);
	Vector y = VectorXd::Zero(m);
	Vector x00 = x0;

#ifdef DIAGONAL_NORMALIZATION
	for (int i = 0; i < dim; i++) {
		double diag_ele = A(i,i);
		b(i) /= diag_ele;
		for (int j = 0; j < dim; j++) {
			A(i, j) /= diag_ele;
		}
	}
#endif // DIAGONAL_NORMALIZATION


	//std::cout << "GMRES Started." << endl;

	while(loop_count++ < kmax){
		//std::cout << "<Loop started:" << setw(3) << loop_count << ">" << endl;

		r0 = b - A * x00;
		if (r0.norm() < ZERO_CHECK) {
			return x00;
		}
		v.col(0) = r0.normalized();
		g(0) = r0.norm();


#ifdef DEBUG_GMRES2
		std::cout << "\nLoop:" << loop_count << endl;
		std::cout << "x00" << endl << x00 << endl;
#endif // DEBUG_GMRES2

		for (int j = 0; j < m; j++) {
			//Arnoldi Process
			Av = A * v.col(j);

			temp_norm = 0;
			temp = Av - v.col(j);
			for (int i = 0; i < dim; i++) temp_norm += temp(i)*temp(i);
			temp_norm = sqrt(temp_norm);
			if(temp_norm < ZERO_CHECK){
				Av_same_as_v = 1;
				break;
			}

			for (int i = 0; i <= j; i++) h(i, j) = Av.dot(v.col(i));
			temp = Av;
			for (int i = 0; i <= j; i++){
				temp -= h(i, j)*v.col(i);
			}

			temp_norm = 0;
			for (int i = 0; i < dim; i++) temp_norm += temp(i)*temp(i);
			temp_norm = sqrt(temp_norm);
			h(j + 1, j) = temp_norm;
			v.col(j + 1) = temp/temp_norm;
		

			//Givens Rotation
			r(0, j) = h(0, j);
			for (int i = 0; i <= j - 1; i++) {
				double d_temp = c(i) * r(i, j) + s(i) * h(i + 1, j);
				r(i + 1, j) = -s(i) * r(i, j) + c(i) * h(i + 1, j);
				r(i, j) = d_temp;
			}
			double delta = sqrt(r(j, j)*r(j, j) + h(j + 1, j)*h(j + 1, j));
			c(j) = r(j, j) / delta;
			s(j) = h(j + 1, j) / delta;
			g(j + 1) = -s(j) * g(j);
			g(j) = c(j)*g(j);
			r(j, j) = c(j)*r(j, j) + s(j)*h(j + 1, j);
		}

		if (Av_same_as_v == 1) {
			Av_same_as_v = 0;
			x00.setRandom();
			std::cout << "GMRES Chages a Initial Vector." << endl;
			loop_count--;
			continue;
		}


		//Back substituition to get y
		double sum_temp = 0;
		for (int i = m-1; i >= 0; i--) {
			for (int j = i + 1; j < m; j++) {
				sum_temp += r(i, j) * y(j);
			}
			y(i) = (g(i) - sum_temp) / r(i, i);
			sum_temp = 0;
		}

		//Calculate approximate solution xm
		xm = x00 + v.block(0,0,dim,m) * y;

#ifdef DEBUG_GMRES1
		std::cout << "Debug" << fixed << endl;
		std::cout << "xm" << endl;
		for (int i = 0; i < dim; i++) {
			std::cout << setprecision(3) << setw(6) << xm(i);
		}std::cout << endl;
		std::cout << "v" << endl;
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < m + 1; j++) {
				std::cout << setprecision(3) << setw(6) << v(i, j);
			}
			std::cout << endl;
		}
		std::cout << "h" << endl;
		for (int i = 0; i < m + 1; i++) {
			for (int j = 0; j < m; j++) {
				std::cout << setprecision(3) << setw(6) << h(i, j);
			}
			std::cout << endl;
		}
		std::cout << "y" << endl;
		for (int i = 0; i < m; i++) {
			std::cout << setprecision(3) << setw(6) << y(i);
		}std::cout << endl;
		std::cout << "r" << endl;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				std::cout << setprecision(3) << setw(6) << r(i, j);
			}
			std::cout << endl;
		}
		std::cout << "g" << endl;
		for (int i = 0; i < m + 1; i++) {
			std::cout << setprecision(3) << setw(6) << g(i);
		}std::cout << endl;
		std::cout << "c" << endl;
		for (int i = 0; i < m; i++) {
			std::cout << setprecision(3) << setw(6) << c(i);
		}std::cout << endl;
		std::cout << "s" << endl;
		for (int i = 0; i < m; i++) {
			std::cout << setprecision(3) << setw(6) << s(i);
		}std::cout << endl;
		std::cout << endl;
#endif // DEGUB_GMRES1

		//Converged then return xm.
		if ((b - A * xm).norm() < tol) {
			return xm;
		} else {
			x00 = xm;
		}
		
	}
	
	std::cout << "Equation didn't converged." << endl;
	exit(1);
	
	return xm;
}

template<class Matrix, class Vector>
Vector GMRES(Matrix &A, Vector &b, Vector &x0, double tol) {
	Vector x = GMRES(A, b, x0, tol, KMAX_FOR_GMRES, M_FOR_GMRES);
	return x;
}