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

using namespace std;

namespace pair_potential {
	class Potential {
	public:
		double energy (const double r) const {return energy_(r) + Ushift_;}
		double force (const double r) const {return force_(r);}
		void set_cutoff (double rcut, double Ushift);
		
	private:
		virtual double energy_ (double r) = 0;			//!< Return U(r)
		virtual vector<double> force_ (double r) = 0;	//!< Return \f[ F(x_i) = -\frac{\del U}{\del r}\frac{\del r}{\del x_i} \f] for each cartesian direction x_i
		double rcut_;									//!< Cutoff radius for interactions
		double Ushift_;									//!< Constant shift to apply to energy function for r < r_{cut}
	};
	
	/*!
	 Set the cutoff values for the potential
	 \param [in] rcut Cutoff radius
	 \param [in] Ushift Energy shift to apply at r \ge r_{cut}
	 */
	void Potential::set_cutoff (double rcut, double Ushift) {
		rcut_ = rcut;
		Ushift_ = Ushift;
	}
	
	//! Shifted Lennard-Jones interaction
	/*!
	 This is the same as standard LJ if \Delta = 0.  This potentially is generally useful for systems with large size asymmetries.
	 \f{eqnarray*}{
	 U(r) &=& 4\epsilon\left(\left(\frac{\sigma}{r-\Delta}\right)^{12} - \left(\frac{\sigma}{r-\Delta}\right)^{6}\right) + U_{shift} & r < r_{cut} \\
	 &=& 0 & r \ge r_{cut}
	 \f}
	 */
	class slj : public Potential {
	public:
		void set_coeff (double eps, double sig, double delta, double rcut);
		
	private:
		double epsilon_;
		double sigma_;
		double delta_;
	};	
	

	/*!
	 \param [in] eps \f[ \epsilon \f] in LJ Equation.
	 \param [in] sig \f[ \sigma \f] in LJ Equation.
	 \param [in] delta \f[ \Delta \f] in LJ Equation.
	 \param [in] rcut \f[ \r_{cut} \f] in LJ Equation.
	 */
	void slj::set_coeff (double eps, double sig, double delta) {
		epsilon_ = eps;
		sigma_ = sig;
		delta_ = delta;
	}
}

#endif