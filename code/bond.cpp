/**
MD Bond Information
\author Nathan A. Mahynski
**/

#include "bond.h"

using namespace bond;
 // George and Arun .. its in the header. not sure if that means it will be autoinitialized? 11/30/12 11:36PM
//char err_msg[ERR_FLAG_SIZE]; //!< Error message buffer commonly used in routines in this namespace
//const double LJ_CUT = pow(2.0, 1/6.0);


// There a number of "standard" bond types, the two that are usually used are FENE and Harmonic, see below

// Fene bond
void Fene::set_coeff (const vector <double> coeffs) {
	assert (coeffs.size() == 5);
	epsilon_ = coeffs[0];
	sigma_ = coeffs[1];
	r0_ = coeffs[2];
	k_ = coeffs[3];
	delta_ = coeffs[4];
}

inline double Fene::energy (const double r2) {
	double r = sqrt(r2), ener = 0.0, b = r-delta_, a = b/r0_;
	if (b < LJ_CUT*sigma_) {
		double x = (sigma_/b), x2 = x*x, x6 = x2*x2*x2, x12 = x6*x6;
		ener += 4.0*epsilon_*(x12-x6)+epsilon_;
	}
	return ener - 0.5*k_*r0_*r0_*log(1-a*a);
}

inline vector <double> Fene::force (const double r2, const double *xyz) {
	vector <double> bond_force(3,0.0);
	double r = sqrt(r2), a = r - delta_;
	for (int i = 0; i < 3; ++i) {
		bond_force[i] = k_*a/(r-a*a*r/(r0_*r0_))*xyz[i];
	}
	if (a < LJ_CUT*sigma_) {
		double x = (sigma_/a), x2 = x*x, x6 = x2*x2*x2, x12 = x6*x6;
		for (int i = 0; i < 3; ++i) {
			bond_force[i] += 24.0*epsilon_/(r*a)*(2.0*x12-x6)*xyz[i];
		}
	}
	return bond_force;
}

// Harmonic Bond
inline double Harmonic::energy (const double r2) {
	double r = sqrt(r2);
	return 0.5*k_*(r-r0_)*(r-r0_);
}

inline vector <double> Harmonic::force (const double r2, const double *xyz) {
	vector <double> bond_force(3,0.0);
	double ir = 1.0/sqrt(r2);
	for (int i = 0; i < 3; ++i) {
		bond_force[i] = -k_*xyz[i]*(1-r0_*ir);
	}
	return bond_force;
}

void Harmonic::set_coeff (const vector <double> coeffs) {
	assert (coeffs.size() == 2);
	r0_ = coeffs[0];
	k_ = coeffs[1];
}
