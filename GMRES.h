#pragma once
#include<Eigen/Dense>
#include<iostream>
using namespace Eigen;
using namespace std;

//#define DEBUG_GMRES

// GMRE5
template<class Matrix, class Vector>
Vector GMRES(Matrix &A, Vector &b, Vector &x0, double tol, int kmax, int m) {
	if ((A.cols() != A.rows()) || (b.size() != x0.size()) || (A.cols() != x0.size())) {
		std::cout << "Matrix, vector size error." << endl;
		exit(EXIT_FAILURE);	
	}

	//Get dimension of matrices.
	int dim = x0.size();
	int loop_count = 0;
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

	while(loop_count++ < kmax){
		std::cout << "<Loop started:" << setw(3) << loop_count << ">" << endl;

		r0 = b - A * x00;
		v.col(0) = r0.normalized();
		g(0) = r0.norm();
		
		for (int j = 0; j < m; j++) {
			//Arnoldi Process
			Av = A * v.col(j);
			for (int i = 0; i <= j; i++) h(i, j) = Av.dot(v.col(i));
			temp = Av;
			for (int i = 0; i <= j; i++){
				temp -= h(i, j)*v.col(i);
			}
			h(j + 1, j) = temp.norm();
			v.col(j + 1) = temp.normalized();
			
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

#ifdef DEBUG_GMRES
		std::cout << "Debug" << fixed << endl;
		std::cout << "xm" << endl;
		for (int i = 0; i < dim; i++) {
			std::cout << setprecision(4) << setw(8) << xm(i);
		}std::cout << endl;
		std::cout << "v" << endl;
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < m + 1; j++) {
				std::cout << setprecision(4) << setw(8) << v(i, j);
			}
			std::cout << endl;
		}
		std::cout << "h" << endl;
		for (int i = 0; i < m + 1; i++) {
			for (int j = 0; j < m; j++) {
				std::cout << setprecision(4) << setw(8) << h(i, j);
			}
			std::cout << endl;
		}
		std::cout << "y" << endl;
		for (int i = 0; i < m; i++) {
			std::cout << setprecision(4) << setw(8) << y(i);
		}std::cout << endl;
		std::cout << "r" << endl;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < m; j++) {
				std::cout << setprecision(4) << setw(8) << r(i, j);
			}
			std::cout << endl;
		}
		std::cout << "g" << endl;
		for (int i = 0; i < m + 1; i++) {
			std::cout << setprecision(4) << setw(8) << g(i);
		}std::cout << endl;
		std::cout << "c" << endl;
		for (int i = 0; i < m; i++) {
			std::cout << setprecision(4) << setw(8) << c(i);
		}std::cout << endl;
		std::cout << "s" << endl;
		for (int i = 0; i < m; i++) {
			std::cout << setprecision(4) << setw(8) << s(i);
		}std::cout << endl;
		std::cout << endl;
#endif // DEGUB_GMRES
		

		//Converged then return xm.
		if ((b - A * xm).norm() < tol) {
			return xm;
		} else {
			x00 = xm;
		}
	}

	std::cout << "Ax=b didn't converged." << endl;
	return xm;
}