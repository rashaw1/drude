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

	// Finally calculate J_R(C)
	JR = 0.0;
	// Calculate the R.C.R part
	for (int i = 0; i < N-1; i++){
		for (int j = 0; j < N-1; j++)
			JR += C(i, j) * ( R(i, 0)*R(j, 0) + R(i, 1)*R(j, 1) + R(i, 2)*R(j, 2) );
	}
	// Exponentiate
	JR = std::exp(-1.0*JR);
}

// Copy constructor
BasisFunction::BasisFunction(const BasisFunction& other)
{
	A = other.A;
	B = other.B;
	C = other.C;
	k = other.k;
	JR = other.JR;
}					 
