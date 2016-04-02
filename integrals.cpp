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
double overlap(int N, Eigen::MatrixXd& A_dp_inv, Eigen::MatrixXd& B_dp, Eigen::MatrixXd& C_dp,
			   Eigen::MatrixXd& R, double detA_dp, double norm1, double norm2, double JR1, double JR2)
{
	// Calculate the K factor
	// K = exp[R.(1/4 B. A''^-1.B - C'').R ] 
	double K = 0.0;
	Eigen::MatrixXd tempMat = 0.25* B_dp * A_dp_inv * B_dp.transpose() - C_dp;
	for (int i = 0; i < N-1; i++){
		for (int j = 0; j < N-1; j++){
			double temp = R(i, 0) * R(j, 0);
			temp += R(i, 1) * R(j, 1);
			temp += R(i, 2) * R(j, 2);
			K += tempMat(i, j) * temp;
		}
	}
	K = K - 0.5*( JR1 + JR2 );
	K = std::exp(K);

	// Calculate the overlap integral
	double S =std::pow(M_PI, N)/detA_dp;
	
	S =  K * std::pow(S, 1.5) * norm1 * norm2;

	if (S < 1e-14) S = 0.0;
		
	return S;
}

// Kinetic energy integral
double kinetic(int N, Eigen::MatrixXd& A, Eigen::MatrixXd& A_p,
			   Eigen::MatrixXd& A_dp, Eigen::MatrixXd& A_dp_inv,
			   Eigen::MatrixXd& k, Eigen::MatrixXd& k_p,
			   Eigen::MatrixXd& k_dp, double S, double mu)
{
	// Calculate the T factor
	// Start with A''^{-1} A^2 A''^{-1}
	Eigen::MatrixXd tempMat = A_dp_inv*A_p*A*A_dp_inv;
	
	// Determine k''. A''^{-1} A' A A''^{-1} . k''
	double T = 0.0;
  	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			double temp = k_dp(i, 0) * k_dp(j, 0) + k_dp(i, 1) * k_dp(j, 1);
			temp += k_dp(i, 2) * k_dp(j, 2);
			T += tempMat(i, j) * temp;
		}
	}

	// Then k . k'
	for (int i = 0; i < N; i++){
			T += k(i, 0) * k_p(i, 0);
			T += k(i, 1) * k_p(i, 1);
			T += k(i, 2) * k_p(i, 2);
	}

	// And -k'.A A''^{-1}.k''
	tempMat = A*A_dp_inv;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			double temp = k_p(i, 0) * k_dp(j, 0);
			temp += k_p(i, 1) * k_dp(j, 1) + k_p(i, 2) * k_dp(j, 2);
			T -= tempMat(i, j) * temp;
		}
	}
	
	// And -k''.A''^-1 A'.k
	tempMat = A_dp_inv * A_p;
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			double temp = k_dp(i, 0) * k(j, 0);
			temp += k_dp(i, 1) * k(j, 1) + k_dp(i, 2) * k(j, 2);
			T -= tempMat(i, j) * temp;
		}
	}

	
	// This just leaves 6Tr{A A''^-1 A'} 
	tempMat = A * A_dp_inv * A_p;
	T += 6.0*tempMat.trace(); 

	T = T*S/(2.0*mu);
	
	return T;
}

// Drude potential integral
double drudePotential(Eigen::MatrixXd& A_dp_inv, Eigen::MatrixXd& B_dp,
					  Eigen::MatrixXd& R, double S, double mu, double omega, int N)
{
   double VD =  1.5 * A_dp_inv.trace();
   Eigen::MatrixXd tempMat = 0.25*B_dp*A_dp_inv*A_dp_inv*B_dp.transpose();

   for (int i = 0; i < N-1; i++){
	   for (int j = 0; j < N-1; j++){
		   double temp =R(i, 0)*R(j, 0) + R(i, 1)*R(j, 1);
		   temp += R(i, 2) * R(j,2);
		   VD += tempMat(i,j)*temp;
	   }
   }

   return 0.5*mu*omega*omega*VD*S;
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
						  Eigen::MatrixXd& R, double mu, double omega, double q,
						  Eigen::MatrixXd& Smat, int s_i, int s_j)
{

	// Construct the double-primed matrices
	Eigen::MatrixXd A_dp, B_dp, C_dp, k_dp;
	A_dp = phi1.getA() + phi2.getA();
	B_dp = phi1.getB() + phi2.getB();
	C_dp = phi1.getC() + phi2.getC();
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
	
	double S = overlap(N, A_dp_inv, B_dp, C_dp, R, detA_dp, phi1.getNorm(), phi2.getNorm(),
					   phi1.getJR(), phi2.getJR());
	Smat(s_i, s_j) = S;
	Smat(s_j, s_i) = S;
	
	// Kinetic energy integral
	double T = kinetic(N, phi1.getA(), phi2.getA(), A_dp, A_dp_inv, phi1.getk(), phi2.getk(),  k_dp, S, mu);

	// Drude potential energy integral
	double VD = drudePotential(A_dp_inv, B_dp, R, S, mu, omega, N);

	// Now calculate the coulombic potential integrals
	double VC = 0.0;
	double Rij = 0.0;
	std::vector<double> Rij_vec(3);
	
	for (int i = 0; i < N-1; i++){
		Rij_vec[0] = R(i, 0); Rij_vec[1] = R(i, 1); Rij_vec[2] = R(i, 2);
		Rij = R(i, 0)*R(i, 0) + R(i, 1)*R(i, 1) + R(i, 2)*R(i, 2);
		Rij = std::sqrt(Rij);

		// Get gi
		double gi =1.0/A_dp_inv(i, i);
		
		for (int j = i+1; j < N; j++){
			// Add the Rij term
			VC += S/Rij;

			std::vector<double> Pa_vec(3), Pb_vec(3), Pc_vec(3);
			double Pa = 0, Pb = 0, Pc = 0;
			for (int m = 0; m < 3; m++){
				Pa_vec[m] = Rij_vec[m];
				Pb_vec[m] = Rij_vec[m];
				Pc_vec[m] = Rij_vec[m];
				for (int n = 0; n < N; n++){
					Pa_vec[m] += 0.5*A_dp_inv(i, n) * k_dp(n, m);
					Pb_vec[m] -= 0.5*A_dp_inv(j, n) * k_dp(n, m);
					Pc_vec[m] -= 0.5*(A_dp_inv(j, n) - A_dp_inv(i, n))*k_dp(n, m);
				}
				Pa += Pa_vec[m]*Pa_vec[m];
				Pb += Pb_vec[m]*Pb_vec[m];
				Pc += Pc_vec[m]*Pc_vec[m];
			}
			Pa = std::sqrt(Pa);
			Pb = std::sqrt(Pb);
			Pc = std::sqrt(Pc);
			
			// Get  gj, and gij
			double gj = 1.0/A_dp_inv(j, j);
			double gij = A_dp_inv(i, i) + A_dp_inv(j, j) - 2.0*A_dp_inv(i, j);
			gij = 1.0/gij;

			// Add the aij, bji, cij terms
			VC += coulombPotential(Pc, gij, S);
			VC -= coulombPotential(Pa, gi, S);
			VC -= coulombPotential(Pb, gj, S);

			// Update Rij
			if (j != N-1) { 
				Rij_vec[0] += R(j, 0); Rij_vec[1] += R(j, 1); Rij_vec[2] += R(j, 2);
				double temp = R(j, 0)*R(j, 0) +R(j, 1)*R(j, 1) + R(j, 2)*R(j, 2);
				Rij += std::sqrt(temp);
			}

		}
	}
	VC = q*q*VC;

	double h_el = T + VD + VC;
	//	std::cout << S << " " << T << " " << VD << " " << VC << "\n";
	
	if (std::isnan(h_el) || std::isinf(h_el)) {  h_el = 0.0; std::cout << "Dodgy hamiltonian element.\n"; }
	
	return h_el;
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
	Eigen::MatrixXd S(nbfs, nbfs);
	for (int i = 0; i < nbfs; i++){
		for (int j = i; j < nbfs; j++){
			H(i, j) = hamiltonianElement(N, bfs[i], bfs[j], R, mu, omega, q, S, i, j);
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
	Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver(H, S, Eigen::EigenvaluesOnly);
	return solver.eigenvalues();
}
		
