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
	
	class Bond {
	public:
		virtual double energy (double r) const = 0;		//!< Pure virtual requirement of an energy function for a bond class 
		virtual double force (double r) const = 0;		//!< Pure virtual requirement of a force function for a bond class
		
	private:
		int bond_index_;								//!< Internal index associated with bond name_ assigned by the user
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
	 U_{WCA} &=& 4\epsilon \left( \left( \frac{\sigma}{r-\Delta} \right)^{12} - \left( \frac{\sigma}{r-\Delta} \right)^6 \right) & r < 2^{1/6}\sigma+\Delta \\
	 &=& 0 &r \ge 2^{1/6}\sigma+\Delta
	 \f}
	 */
	class Fene : public Bond {
	public:
		double energy (double r);
		double force (double r);
		void set_coeff (double eps, double sig, double r0, double k); //!< Set the coefficients for the bond potential
		
	private:
		double epsilon_;							//!< Repulsive force strength (in energy units)
		double sigma_;								//!< Repulsive force interaction distance (in distance units)
		double r0_;									//!< Size parameter (in distance units)
		double k_;									//!< Attractive force strength (in units of energy/distance^2)
		double delta_;								//!< Shift to account for molecules of different diameters, \f[ \Delta = (d_i + d_j)/2 - 1 \f], this is automatically calculated
	};
	
	void Fene::set_coeff (double eps, double sig, double r0, double k) {
		epsilon_ = eps;
		sigma_ = sig;
		r0_ = r0;
		k_ = k;
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
		double energy (double r);
		double force (double r);
		void set_coeff (double r0, double k);		//!< Set the coefficients for the bond potential
		
	private:
		double r0_;									//!< Minimum of energy, equilibrium location
		double k_;									//!< Bond force strength (in units of energy/distance^2)
	};
	
	void Harmonic::set_coeff (double r0, double k) {
		r0_ = r0;
		k_ = k;
	}
	
}

#endif