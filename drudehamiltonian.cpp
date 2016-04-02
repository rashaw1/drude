/*************************************************************
 *
 * PURPOSE: To construct and diagonalise the Hamiltonian for 
 *          a system of N drude oscillators.
 *
 * AUTHOR: Robert shaw
 *************************************************************/

#include <iostream>
#include <iomanip>
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
	double base2 = 2.0;
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
	std::cout << "R:\n";
	std::cin >> R(0, 0);
	std::cout << "Zeta:\n";
	std::cin >> zetaval;
	
	std::cout << "Min. exp:\n";
	std::cin >> minexp;
	std::cout << "No. exps:\n";
	std::cin >> nexps;
	std::cout << "Exponent base:\n";
	std::cin >> base;
	std::cout << "Zeta base:\n";
	std::cin >> base2;
	
    // Initialise array of basis functions
	std::vector<BasisFunction> bfs;
	
	// Make parameter matrices/vectors
	std::vector<double> zeta(N, zetaval);
	Eigen::MatrixXd alpha(N, N);
	Eigen::MatrixXd beta(N, N);
	Eigen::MatrixXd gamma(N, N);

	int nbfs = 0;
	
	// Make basis functions
	for (int i = 0; i < nexps+1; i++){

		alpha(0, 1) = std::pow(base, minexp + i);
		if (i == nexps)  alpha(0, 1) = 0.0; 
		
		for (int j = 0; j < nexps+ 1; j++){

			beta(1, 0) = std::pow(base, minexp + j);
			if (j == nexps) beta(1, 0) = 0.0;
			
			for(int k = 0; k < nexps+1; k++){

				gamma(0, 1) = std::pow(base, minexp + k);
				if (k == nexps) gamma(0, 1) = 0.0;

				for (int m = 0; m < 1; m++){
					zeta[0] = zetaval*std::pow(base2, m);

					for (int n = 0; n < 1; n++){
						zeta[1] = zetaval*std::pow(base2, n);
				
						BasisFunction bftemp(N, zeta, alpha, beta, gamma, R);
						if(!bftemp.istoosmall()){
							bfs.push_back(bftemp);
							nbfs++;
						}
					}
				}
			}
		}
	}
	std::cout << "Total number of basis functions = " << nbfs << "\n";
	
	// Form and diagonalise hamiltonian matrix
	Eigen::MatrixXd D = hamiltonian(N, nbfs, bfs, R, mu, omega, q); 

	std::cout << "Hamiltonian formed and solved.\n";

	int i = 0;
	while ( D(i) < 0.1 ) i++;
  	std::cout << "Lowest eigenvalue = " << std::setprecision(15) << D(i) << "\n" ;

	return 0;
}

