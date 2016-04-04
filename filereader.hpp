/*************************************************************
 *
 * PURPOSE: To handle input and basis file reading. 
 *
 * AUTHOR: Robert Shaw
 *
 *************************************************************/

#ifndef FILEREADERHDRDEF
#define FILEREADERHDRDEF

#include <fstream>
#include <vector>
#include <Eigen/Dense>

class BasisFunction;

// Read in the number of oscillators, and their parameters, mu, omega, q
int readParams(std::ifstream& input, std::vector<double>& mu,
				std::vector<double>& omega, std::vector<double>& q);

// Read in the geometry matrix R
void readGeom(std::ifstream& input, Eigen::MatrixXd& R, int N);

// Read in the basis set
int readBasis(std::ifstream& basis, std::vector<BasisFunction>& bfs, int N,
			  std::vector<double>& zeta, Eigen::MatrixXd& R);


#endif
