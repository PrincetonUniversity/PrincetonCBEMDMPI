/**
 MD Bond Information
 \author Nathan A. Mahynski
 **/

#ifndef BOND_H_
#define BOND_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include "misc.h"

using namespace std;
using namespace misc;

//! Namespace containing pertinent Bond information
namespace bond {
	char err_msg[ERR_FLAG_SIZE]; //!< Error message buffer commonly used in routines in this namespace
	const double LJ_CUT = pow(2.0, 1/6.0);
	
	class Bond {
	public:
		virtual inline double energy (const double r2) const = 0;								//!< Pure virtual requirement of an energy function for a bond class 
		virtual inline vector <double> force (const double r2, const double *xyz) const = 0;	//!< Pure virtual requirement of a force function for a bond class
		virtual inline void set_coeff (const vector <double> coeffs) = 0;						//!< Set the coefficients for a bond potential
		int bond_index;																			//!< Internal index associated with bond name assigned by the user
	};
	
	// There a number of "standard" bond types, the two that are usually used are FENE and Harmonic, see below
	
	//! Finitely Extensible Non-linear Elastic Bond
	/*!
	 The Fene bond potential is given by:
	 \f[
	 U(r) = -\frac{1}{2}kr_{0}^2\text{ln}\left(1-\left(\frac{r-\Delta}{r_0}\right)^2\right) + U_{WCA}
	 \f]
	 Where the short range repulsion is provided by the WCA potential:
	 \f{eqnarray*}{
	 U_{WCA} &=& 4\epsilon \left( \left( \frac{\sigma}{r-\Delta} \right)^{12} - \left( \frac{\sigma}{r-\Delta} \right)^6 \right) + \epsilon & r < 2^{1/6}\sigma+\Delta \\
	 &=& 0 &r \ge 2^{1/6}\sigma+\Delta
	 \f}
	 */
	class Fene : public Bond {
	public:
		inline double energy (const double r2);
		inline vector <double> force (const double r2, const double *xyz);
		inline void set_coeff (const vector <double> coeffs);
		
	private:
		double epsilon_;							//!< Repulsive force strength (in energy units)
		double sigma_;								//!< Repulsive force interaction distance (in distance units)
		double r0_;									//!< Size parameter (in distance units)
		double k_;									//!< Attractive force strength (in units of energy/distance^2)
		double delta_;								//!< Shift to account for molecules of different diameters, \f[ \Delta = (d_i + d_j)/2 - 1 \f], this is automatically calculated
	};
	
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
		double a = r - delta_;
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
	
	//! Harmonic Bond
	/*!
	 The Harmonic bond potential is given by:
	 \f[
	 U(r) = \frac{1}{2}k\left(r - r_{0}\right)^2
	 \f]
	 */
	class Harmonic : public Bond {
	public:
		inline double energy (const double r2);
		inline vector <double> force (const double r2, const double *xyz);
		inline void set_coeff (const vector <double> coeffs);		
		
	private:
		double r0_;									//!< Minimum of energy, equilibrium location
		double k_;									//!< Bond force strength (in units of energy/distance^2)
	};
	
	inline double Harmonic::energy (const double r2) {
		double r = sqrt(r2);
		return 0.5*k_*(r-r0_)*(r-r0_);
	}
	
	inline vector <double> Harmonic::force (const double r2, const double *xyz); {
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
}

#endif