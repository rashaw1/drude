/*************************************************************
 *
 * PURPOSE: To implement the file reading capabilities
 *          declared in filereader.hpp
 *
 * AUTHOR: Robert Shaw
 *
 *************************************************************/

#include "filereader.hpp"
#include "basisfunction.hpp"
#include <iomanip>
#include <algorithm>
#include <iostream>

int readParams(std::ifstream& input, std::vector<double>& mu,
				std::vector<double>& omega, std::vector<double>& q)
{
	std::string line, token;
	std::size_t pos;
	int N = 0;
	// Get line and parse
	while ( std::getline(input, line) ) {
		pos = line.find(',');
		if ( pos != std::string::npos ) {

			// Tokenise and remove spaces
			token = line.substr(0, pos);
			token.erase(std::remove(token.begin(), token.end(), ' '), token.end());
			// Make lowercase
			std::transform(token.begin(), token.end(), token.begin(), ::tolower);
			
			// Find the token
			if ( token == "n" ) { // Number of oscillators
				N = std::stoi( line.substr(pos+1, line.length() ));
			} else if (token == "mu" ) {
				
				line = line.substr( pos+1, line.length() );
				double mu_temp;

				std::size_t pos2 = line.find(',');
				while ( pos2 != std::string::npos ){
					mu_temp = std::stod( line.substr(0, pos2) );
					mu.push_back(mu_temp);
				   
					line = line.substr( pos2+1, line.length() );
					pos2 = line.find(',');
				}
				mu_temp = std::stod( line );
				mu.push_back(mu_temp);
			} else if (token == "omega") {
				line = line.substr( pos+1, line.length() );
				double omega_temp;

				std::size_t pos2 = line.find(',');
				while ( pos2 != std::string::npos ){
					omega_temp = std::stod( line.substr(0, pos2) );
					omega.push_back(omega_temp);
					pos = pos2;
					line = line.substr( pos2 + 1, line.length() );
					pos2 = line.find(',');
				}
				omega_temp= std::stod( line );
				omega.push_back(omega_temp);
			} else if (token == "q") {
				line = line.substr( pos+1, line.length() );
				double q_temp;

				std::size_t pos2 = line.find(',');
				while ( pos2 != std::string::npos ){
					q_temp = std::stod( line.substr(0, pos2) );
					q.push_back(q_temp);
					pos = pos2;
					line = line.substr( pos2 + 1, line.length() );
					pos2 = line.find(',');
				}
				q_temp= std::stod( line );
				q.push_back(q_temp);
			}
		}
	}

	// Check vector sizes
	if ( mu.size() != N ) 
		std::cerr << "Wrong number of mu parameters.\n";
	else if ( omega.size() != N )
		std::cerr << "Wrong number of omega parameters.\n";
	else if ( q.size() != N )
		std::cerr << "Wrong number of q parameters.\n";

	return N;
}

// Read the geometry
void readGeom(std::ifstream& input, Eigen::MatrixXd& R, int N)
{
	std::string line;
	bool found = false;
   
	while ( std::getline(input, line) && !found){

		std::size_t pos = line.find(',');
		if ( line.substr(0, pos) == "geom" ) found = true;
			
	}

	if (!found) std::cerr << "Geometry not found.\n";
	else {
		// Read the geometry
		for (int i = 0; i < N-1; i++){
			std::size_t pos = line.find(',');
			R(i, 0) = std::stod( line.substr(0, pos) );
			line = line.substr(pos+1, line.length() );
			pos = line.find(',');
			R(i, 1) = std::stod( line.substr(0, pos) );
			R(i, 2) = std::stod( line.substr(pos+1, line.length()) );
			std::getline(input, line);
		}
	}
}

// Read in the basis functions
int readBasis(std::ifstream& basis, std::vector<BasisFunction>& bfs, int N,
			  std::vector<double>& zeta, Eigen::MatrixXd& R)
{
	int nbfs = 0;
	std::string line;
	std::size_t pos;
	Eigen::MatrixXd alpha(N, N);
	Eigen::MatrixXd beta(N, N);
	Eigen::MatrixXd gamma(N, N);

	while ( std::getline(basis, line) ){
		// Remove spaces
		line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
		// Make lowercase
		std::transform(line.begin(), line.end(), line.begin(), ::tolower);

		if ( line == "basisfunction" ){
			bool endfunc = false;
			while (!endfunc ){
				std::getline(basis, line);
				// Remove spaces
				line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
				// Make lowercase
				std::transform(line.begin(), line.end(), line.begin(), ::tolower);
				if ( line == "alpha," ) {
					for (int i = 0; i < N; i++){
						std::getline(basis, line);
						for (int j = 0; j < N; j++){
							pos = line.find(',');
							alpha(i, j) = std::stod( line.substr(0, pos) );
							line = line.substr( pos+1, line.length() );
						}
					}
				} else if (line == "beta,") {
					for (int i = 0; i < N; i++){
						std::getline(basis, line);
						for (int j = 0; j < N; j++){
							pos = line.find(',');
							beta(i, j) = std::stod( line.substr(0, pos) );
							line = line.substr( pos+1, line.length() );
						}
					}
				} else if (line == "gamma,") {
					for (int i = 0; i < N; i++){
						std::getline(basis, line);
						for (int j = 0; j < N; j++){
							pos = line.find(',');
							gamma(i, j) = std::stod( line.substr(0, pos) );
							line = line.substr( pos+1, line.length() );
						}
					}
				} else if (line == "endfunction") {
					BasisFunction bftemp(N, zeta, alpha, beta, gamma, R);
					bfs.push_back(bftemp);
					nbfs++;
					endfunc = true;
				}
			}
		}
	}
	return nbfs;
}
