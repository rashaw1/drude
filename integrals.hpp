/*************************************************************
 *
 * PURPOSE: To declare functions for calculating the matrix
 *          elements between two BasisFunctions.
 *
 * AUTHOR: Robert Shaw
 *
 *************************************************************/

#ifndef INTHDRDEF
#define INTHDRDEF

// Forward declare BasisFunction
#include "basisfunction.hpp"
#include <Eigen/Dense> 
#include <Eigen/Eigenvalues>
#include <vector>

// Overlap integral
// Takes the no. of oscillators, the matrices A''^-1 and k'' = k + k' as
// arguments, along with the two factors J_R(C) and J_R(C'),
// and the determinant of A''.
double overlap(int N, Eigen::MatrixXd& A_dp_inv, Eigen::MatrixXd& B_dp, Eigen::MatrixXd& C_dp,
			   Eigen::MatrixXd& R, double detA_dp, double norm1, double norm2, double JR1, double JR2);

// Kinetic energy integral
// Takes that matrices A, A'', its inverse A''^-1, k, and k'' as arguments,
// along with the overlap integral, S, and the parameter matrix L2
double kinetic(int N, Eigen::MatrixXd& A, Eigen::MatrixXd& A_p,
			   Eigen::MatrixXd& A_dp, Eigen::MatrixXd& A_dp_inv,
				   Eigen::MatrixXd& k, Eigen::MatrixXd& k_p,
			   Eigen::MatrixXd& k_dp, double S, Eigen::MatrixXd& L2);

// Drude potential energy integral
// Takes the matrix A''^-1 and the overlap integral S as arguments
// along with the parameters mu and omega
double drudePotential(Eigen::MatrixXd& A_dp_inv, Eigen::MatrixXd& B_dp, Eigen::MatrixXd& R,
					  double S, Eigen::MatrixXd& M, int N);

// Coulomb potential energy integrals
// Takes the quantities Rij, g, and the overlap integral S as arguments
double coulombPotential(double Rij, double g, double S);

// Calculate the Hamiltonian matrix element
// between BasisFunctions phi1 and phi2.
// Takes the number of oscillators, N,
// basis functions, and the matrix R as arguments,
// along with the drude parameters, mu, omega, and q
double hamiltonianElement(int N, BasisFunction& phi1, BasisFunction& phi2,
						  Eigen::MatrixXd& R, Eigen::MatrixXd& L2, Eigen::MatrixXd& M,
						  std::vector<double> q, Eigen::MatrixXd& Smat, int i, int j);

// Construct the Hamiltonian matrix
Eigen::MatrixXd hamiltonian(int N, int nbfs,
							std::vector<BasisFunction> bfs, Eigen::MatrixXd& R,
							std::vector<double> mu, std::vector<double> omega, std::vector<double> q);

#endif
