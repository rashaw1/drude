/*************************************************************
 * 
 * PURPOSE: To implement integrals.hpp
 *
 * AUTHOR: Robert Shaw
 *
 *************************************************************/

#include "integrals.hpp"
#include <iostream>
#include <cmath>

// Overlap integral
double overlap(int N, Eigen::MatrixXd& A_dp_inv, Eigen::MatrixXd& k_dp,
			   double JRC, double JRCP, double detA_dp)
{
	// Calculate the K factor
	// K = exp[1/4 k''. A''^-1 . k''] 
	double K = 0.0;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			double temp = k_dp(i, 0) * k_dp(j, 0);
			temp += k_dp(i, 1) * k_dp(j, 1);
			temp += k_dp(i, 2) * k_dp(j, 2);
			K += A_dp_inv(i, j) * temp;
		}
	}
	K = std::exp(0.25*K);

	// Calculate the overlap integral
	double S = std::pow(M_PI, N)/detA_dp;

	S = JRC * JRCP * K * std::pow(S, 1.5);

	return S;
}

// Kinetic energy integral
double kinetic(int N, Eigen::MatrixXd& A, Eigen::MatrixXd& A_dp,
			   Eigen::MatrixXd& A_dp_inv, Eigen::MatrixXd& k,
			   Eigen::MatrixXd& k_dp, double S, double mu)
{
	// Calculate the T factor
	// Firstly k'' . A''^-1 A^2 A''^-1 . k''
	Eigen::MatrixXd tempMat = A_dp_inv * A * A * A_dp_inv;
	double T = 0.0;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			double temp = k_dp(i, 0) * k_dp(j, 0);
			temp += k_dp(i, 1) * k_dp(j, 1);
			temp += k_dp(i, 2) * k_dp(j, 2);
			T += tempMat(i, j) * temp;
		}
	}

	// Then -2k. A A''^-1 . k''
	tempMat = A * A_dp_inv;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			double temp = k(i, 0) * k_dp(j, 0);
			temp += k(i, 1) * k_dp(j, 1);
			temp += k(i, 2) * k_dp(j, 2);
			T += -2.0 * tempMat(i, j) * temp;
		}
	}

	// Then k^2
	for (int i = 0; i < N; i++){
			T += k(i, 0) * k(i, 0);
			T += k(i, 1) * k(i, 1);
			T += k(i, 2) * k(i, 2);
	}

	// Finally, -6Tr{A}
	T -= 6.0*A.trace();

	// This just leaves 6Tr{A A''^-1 A}
	tempMat = A * A_dp_inv * A;
	T += 6.0*tempMat.trace();

	T = -T*S/(2.0*mu);

	return T;
}

// Drude potential integral
double drudePotential(Eigen::MatrixXd& A_dp_inv, double S,
					  double mu, double omega)
{
   return 0.75*mu*omega*omega * A_dp_inv.trace();
}

// Coulomb potential integral
double coulombPotential(double Rij, double g, double S)
{
	double V = std::erf(std::sqrt(g) * Rij);
	V = V*S/Rij;

	return V;
}

// Get the total hamiltonian matrix element between phi1 and phi2
double hamiltonianElement(int N, BasisFunction& phi1, BasisFunction& phi2,
				   Eigen::MatrixXd& R, double mu, double omega, double q)
{

	// Construct the double-primed matrices
	Eigen::MatrixXd A_dp, k_dp;
	A_dp = phi1.getA() + phi2.getA();
	k_dp = phi1.getk() + phi2.getk();

	// Diagonalise A''
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(A_dp);
	Eigen::MatrixXd U = solver.eigenvectors(); // Eigenvectors
	Eigen::MatrixXd D = solver.eigenvalues().asDiagonal(); // Eigenvalues

	// Compute determinant and inverse of A''
	Eigen::MatrixXd D_inv = D;
	double detA_dp = 1.0;
	for (int i = 0; i < N; i++){
		D_inv(i, i) = 1.0/D(i, i);
		detA_dp *= D(i, i);
	}
	Eigen::MatrixXd A_dp_inv = U * D_inv * U.transpose();

	// Calculate overlap integral
	double S = overlap(N, A_dp_inv, k_dp, phi1.getJR(), phi2.getJR(), detA_dp);
	
	// Kinetic energy integral
	double T = kinetic(N, phi1.getA(), A_dp, A_dp_inv, phi1.getk(), k_dp, S, mu); 

	// Drude potential energy integral
	double VD = drudePotential(A_dp_inv, S, mu, omega);

	// Now calculate the coulombic potential integrals
	double VC = 0.0;
	double Rij = 0.0;
	
	for (int i = 0; i < N-1; i++){
		Rij = R(i, 0)*R(i, 0) + R(i, 1)*R(i, 1) + R(i, 2)*R(i, 2);
		Rij = std::sqrt(Rij);

		// Get gi
		double gi =1.0/A_dp_inv(i, i);

		for (int j = i+1; j < N; j++){
			// Add the Rij term
			VC += S/Rij;

			// Get  gj, and gij
			double gj = 1.0/A_dp_inv(j, j);
			double gij = A_dp_inv(i, i) + A_dp_inv(j, j) - 2.0*A_dp_inv(i, j);
			gij = 1.0/gij;

			// Add the aij, bji, cij terms
			VC += coulombPotential(Rij, gij, S);
			VC -= coulombPotential(Rij, gi, S);
			VC -= coulombPotential(Rij, gj, S);

			// Update Rij
			if (j != N-1) { 
				double temp = R(j, 0)*R(j, 0) +R(j, 1)*R(j, 1) + R(j, 2)*R(j, 2);
				Rij += std::sqrt(temp);
			}

		}
	}
	VC = q*q*VC;

	return (T + VD + VC);															 
}

// Construct the Hamiltonian matrix
Eigen::MatrixXd hamiltonian(int N, int nbfs,
						   std::vector<BasisFunction> bfs, Eigen::MatrixXd& R,
							double mu, double omega, double q)
{
	int nelements = ( nbfs*(nbfs+1) ) / 2;
	int onepercent = nelements/100;

	std::cout << "Progress:\n";
	int ndone = 0;
	
	Eigen::MatrixXd H(nbfs, nbfs);
	for (int i = 0; i < nbfs; i++){
		for (int j = i; j < nbfs; j++){
			H(i, j) = hamiltonianElement(N, bfs[i], bfs[j], R, mu, omega, q);
			H(j, i) = H(i, j);
			ndone++;
			if ( ndone == onepercent ){
				std::cout << "|";
				std::cout.flush();
				ndone = 0;
			}
		}
	}

	std::cout << "\n\nDiagonalising...\n\n";
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(H);
	return solver.eigenvalues().asDiagonal();
}
		
