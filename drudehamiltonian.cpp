/*************************************************************
 *
 * PURPOSE: To construct and diagonalise the Hamiltonian for 
 *          a system of N drude oscillators.
 *
 * AUTHOR: Robert shaw
 *************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include "basisfunction.hpp"
#include "integrals.hpp"

int main(int argc, char* argv[])
{

	// Parameters
	int N = 2;
	double mu = 0.3020;
	double omega = 0.7272;
	double q = 1.3314;
    int minexp = -2; // minimum exponent
	int nexps = 2; // number of exponents
	double base = 2.0;
	double zetaval = 1.0;
	
	// Geometry
	Eigen::MatrixXd R(N-1, 3);
	R(0, 1) = 0.0; R(0, 2) = 0.0;

	// Get input from user
	std::cout << "Mu:\n";
	std::cin >> mu;
	std::cout << "Omega:\n";
	std::cin >> omega;
	std::cout << "q:\n";
	std::cin >> q;
	std::cout << "Enter R:\n";
	std::cin >> R(0, 0);

	std::cout << "Min. exp:\n";
	std::cin >> minexp;
	std::cout << "No. exps:\n";
	std::cin >> nexps;
	std::cout << "Exponent base:\n";
	std::cin >> base;
	
    // Initialise array of basis functions
	std::vector<BasisFunction> bfs;
	
	// Make parameter matrices/vectors
	std::vector<double> zeta(N, zetaval);
	Eigen::MatrixXd alpha(N, N);
	Eigen::MatrixXd beta(N, N);
	Eigen::MatrixXd gamma(N, N);

	// Make basis functions
	for (int i = 0; i < nexps; i++){

		alpha(0, 1) = std::pow(base, minexp + i);

		for (int j = 0; j < nexps; j++){

			beta(1, 0) = std::pow(base, minexp + j);

			for(int k = 0; k < nexps; k++){

				gamma(0, 1) = std::pow(base, minexp + k);

				BasisFunction bftemp(N, zeta, alpha, beta, gamma, R);
				bfs.push_back(bftemp);
			}
		}
	}

	// Form and diagonalise hamiltonian matrix
	Eigen::MatrixXd D = hamiltonian(N, nexps*nexps*nexps, bfs, R, mu, omega, q); 

	std::cout << "Hamiltonian formed and solved.\n";
	std::cout << "Lowest eigenvalue = " << D(0, 0) << "\n"; 
	
	return 0;
}
