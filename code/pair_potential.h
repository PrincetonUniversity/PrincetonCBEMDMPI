/**
 Pair Potentials
 \author Nathan A. Mahynski
 **/

#ifndef PPOT_H_
#define PPOT_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include "misc.h"

using namespace std;

namespace pair_potential {
	class Potential {
	public:
		// If distances were alredy computed, no need to recompute
		void set_cutoff (const double nrcut, const double nUshift);
		virtual void set_coeff (const vector <double> coeffs) = 0;								//!< Set coefficients in pair potential
		virtual inline double energy (const double r2) = 0;										//!< Returns U(r)
		virtual inline vector <double> force (const double r2, const double *xyz) = 0;			//!< Computes \f[ F(x_i) = -\frac{\del U}{\del r}\frac{\del r}{\del x_i} \f] for each cartesian direction x_i
		
		double rcut;									//!< Cutoff radius for interactions
		double Ushift;									//!< Constant shift to apply to energy function for r < r_{cut}
	};
	
	
	//! Shifted Lennard-Jones interaction
	/*!
	 This is the same as standard LJ if \Delta = 0.  This is generally useful for systems with large size asymmetries.
	 \f{eqnarray*}{
	 U(r) &=& 4\epsilon\left(\left(\frac{\sigma}{r-\Delta}\right)^{12} - \left(\frac{\sigma}{r-\Delta}\right)^{6}\right) + U_{shift} & r - \Delta < r_{cut} \\
	 &=& 0 & r - \Delta \ge r_{cut}
	 \f}
	 */
	class slj : public Potential {
	public:
		void set_coeff (const vector <double> coeffs);
		inline double energy (const double r2);
		//inline void force (Atom *a1, Atom *a2, const double r2, const double *xyz);	
		inline vector <double> force (const double r2, const double *xyz);
		
	private:
		double epsilon_;
		double sigma_;
		double delta_;
	};	

}

#endif
