/*************************************************************
 *
 * PURPOSE: To implement the BasisFunction class.
 *
 * AUTHOR: Robert Shaw
 *
 *************************************************************/

#include "basisfunction.hpp"
#include <iostream>
#include <cmath>
#include <Eigen/Eigenvalues>

// Constructor
BasisFunction::BasisFunction(int N, std::vector<double> zeta,
							 Eigen::MatrixXd alpha, Eigen::MatrixXd beta,
							 Eigen::MatrixXd gamma, Eigen::MatrixXd R)
{
	// Initialise the matrices A, B, C, k to the correct sizes
	A.resize(N, N);
	B.resize(N-1, N);
	C.resize(N-1, N-1);
	k.resize(N, 3);

	// Determine the elements of A
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){

			// Diagonal elements
			if ( i == j ) {
				A(i, j) = zeta[i];
				for (int l = i+1; l < N; l++)
					A(i, i) += alpha(i, l) + gamma(i, l);
				for (int l = 0; l < i; l++)
					A(i, i) += beta(i, l) + gamma(l, i);
			}
			// Off diagonal
			else if ( i < j ) {
				A(i, j) = -gamma(i, j);
			} else {
				A(i, j) = -gamma(j, i);
			}

		}
	}
	
    // Now B
	for (int l = 0; l < N-1; l++){

		// i <= l
		for (int i = 0; i < l+1; i++){
			B(l, i) = 0.0;
			for (int j = l + 1; j < N; j++)
				B(l, i) += alpha(i, j) + gamma(i, j);
			B(l, i) *= -2.0;
		}

		// i > l
		for (int i = l+1; i < N; i++){
			B(l, i) = 0.0;
			for (int j = 0; j < l+1; j++)
				B(l, i) += beta(i, j) + gamma(j, i);
			B(l, i) *= 2.0;
		}
	}

	// Finally, C
	for (int m = 0; m < N-1; m++) {
		for (int n = 0; n < N-1; n++){

			int maxmn = (m ? m > n : n);
			int minmn = (m ? m < n : n);

			C(m, n) = 0.0;
			
			for (int i = 0; i < minmn+1; i++){
				for (int j = maxmn+1; j < N; j++){
					C(m, n) += alpha(i, j) + gamma(i, j) + beta(j, i);
				}
			}
		}
	}

	// Then k = B^T R
	k = B.transpose()*R;

	// Calculate norm
	toosmall = false;

	calcNorm(R, N);
	
}

// Copy constructor
BasisFunction::BasisFunction(const BasisFunction& other)
{
	A = other.A;
	B = other.B;
	C = other.C;
	k = other.k;
	JR = other.JR;
	norm = other.norm;
}

void BasisFunction::calcNorm(Eigen::MatrixXd& R, int N)
{
	// Calculate A_inv and detA
    // Diagonalise 2A
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(2.0*A);
	Eigen::MatrixXd U = solver.eigenvectors(); // Eigenvectors
	Eigen::MatrixXd D = solver.eigenvalues().asDiagonal(); // Eigenvalues

	
	// Compute determinant and inverse of 2A
	Eigen::MatrixXd D_inv = D;
	double detA = 1.0;
	for (int i = 0; i < N; i++){
		D_inv(i, i) = 1.0/D(i, i);
		detA *= D(i, i);
	}
	Eigen::MatrixXd A_inv = U * D_inv * U.transpose();

	// Calculate the K factor
	// K = exp[R.(1/4 2B. (2A)^-1 . 2B - 2C).R]
	double K = 0.0;
	Eigen::MatrixXd tempMat = B * A_inv * B.transpose() - 2.0*C;
	for (int i = 0; i < N-1; i++){
		for (int j = 0; j < N-1; j++){
			double temp = R(i, 0) * R(j, 0);
			temp +=  R(i, 1) * R(j, 1);
			temp += R(i, 2) * R(j, 2);
			K += tempMat(i, j) * temp;
		}
	}
	/*	JR = K;
		
	// Calculate the overlap integral
	norm = std::pow(M_PI, N)/detA;
	norm = std::pow(norm, 0.75);
	norm = 1.0/norm;
*/
				norm = 1.0;
				JR = 0.0;
}

double BasisFunction::getNorm()
{
	return norm;
}
