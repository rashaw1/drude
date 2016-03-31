/*************************************************************
 *
 * PURPOSE: To define a class to wrap a basis function for a
 *          system of N drude oscillators. 
 * 
 * AUTHOR: Robert Shaw
 *
 *************************************************************/

#ifndef BFHDRDEF
#define BFHDRDEF

#include <eigen/dense>
#include <vector>

class BasisFunction
{
private:
	Eigen::MatrixXd A, B, C, k;
	double JR;
	double norm;
	bool toosmall;
public:
	// Constructor
	BasisFunction() { }; // Default constructor
	BasisFunction(int N, std::vector<double> zeta, Eigen::MatrixXd alpha,
				  Eigen::MatrixXd beta, Eigen::MatrixXd gamma,
				  Eigen::MatrixXd R);
	BasisFunction(const BasisFunction& other); // Copy constructor

	// Calculate norm
	void calcNorm(Eigen::MatrixXd& R, int N);
	
	// Accessors
	Eigen::MatrixXd& getA() { return A; }
	Eigen::MatrixXd& getB() { return B; }
	Eigen::MatrixXd& getC() { return C; }
	Eigen::MatrixXd& getk() { return k; }
	double getJR() { return JR; }
	double getNorm();
	bool istoosmall() { return toosmall; }
};

#endif 
