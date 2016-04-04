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
#include "filereader.hpp"

int main(int argc, char* argv[])
{

	int program = 0;

	// Get input file from first command line argument
	if(argc < 3){ // No input file given
		std::cerr << "Usage: ./drudeh [input_file] [basis_file]\n";
		program = -1;
	} else {
		std::string ifname = argv[1]; // Input filename
		std::string bfname = argv[2]; // Basis filename
		std::string ofname = ifname; // Output file prefix
		std::size_t pos = ofname.find('.');
		if (pos != std::string::npos) { // Cut off extension
			ofname.erase(pos, ofname.length());
		}
		ofname += ".output";

		// Open the input file
		std::ifstream input(ifname);
		// Check it opened successfully
		if (!input.is_open()){
			std::cerr << "Failed to open input file.\n";
			program = -1;
		} else {

			// Declare and read parameters
			int N;
			std::vector<double> mu, omega, q;
			N = readParams(input, mu, omega, q);

			// Make the zeta vector
			std::vector<double> zeta(N);
			for (int i = 0; i < N; i++) zeta[i] = 0.5*mu[i]*omega[i];

			// Geometry
			Eigen::MatrixXd R(N-1, 3);
			// Rewind input file
			input.clear();
			input.seekg(0, std::ios::beg);
			// Read in geometry
			readGeom(input, R, N);
						
			// Open basis file
			std::ifstream basis(bfname);
			// Check if opened successfully
			if (!basis.is_open()){
				std::cerr << "Failed to open basis file.\n";
				program = -1;
			} else {
				
				// Initialise array of basis functions
				std::vector<BasisFunction> bfs;
								
				// Read in the basis functions
				int nbfs = readBasis(basis, bfs, N, zeta, R);

				// Form and diagonalise hamiltonian matrix
				Eigen::MatrixXd D = hamiltonian(N, nbfs, bfs, R, mu, omega, q);

				// Find lowest non-zero eigenvalue
				int i = 0;
				double lowest_eig = 0.0;
				while ( D(i) < 0.1 ) i++;
				if ( i < nbfs ) 
					lowest_eig = D(i);  
								
				// Open output file
				std::ofstream output(ofname);
				if (!output.is_open()){
					std::cout << "Couldn't open output file.\n";
					std::cout << "Total number of basis functions = " << nbfs << "\n";
					std::cout << "Lowest eigenvalue = " << std::setprecision(15) << lowest_eig << "\n";
				} else {
					output << "Total number of basis functions = " << nbfs << "\n";
					output << "Lowest eigenvalue = "<< std::setprecision(15) << lowest_eig << "\n";
 				}
			}
		}
	}
	return program;
}

