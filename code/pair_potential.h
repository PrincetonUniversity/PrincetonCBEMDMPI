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
		// If distances were alredy computed, no need to recompute
		void set_cutoff (const double nrcut, const double nUshift);
		virtual void set_coeff (const vector <double> coeffs) = 0;						//!< Set coefficients in pair potential
		virtual inline double energy (const double r2);									//!< Return U(r)
		virtual inline void force (Atom *a1, Atom *a2, const double r2, const double *xyz);	//!< Stores \f[ F(x_i) = -\frac{\del U}{\del r}\frac{\del r}{\del x_i} \f] for each cartesian direction x_i on a1 and a2
		
		double rcut;									//!< Cutoff radius for interactions
		double Ushift;									//!< Constant shift to apply to energy function for r < r_{cut}
	};
	
	/*!
	 Set the cutoff values for the potential
	 \param [in] rcut Cutoff radius
	 \param [in] Ushift Energy shift to apply at r \ge r_{cut}
	 */
	void Potential::set_cutoff (const double nrcut, const double nUshift) {
		rcut = rcut;
		Ushift = Ushift;
	}
	
	//! Shifted Lennard-Jones interaction
	/*!
	 This is the same as standard LJ if \Delta = 0.  This is generally useful for systems with large size asymmetries.
	 \f{eqnarray*}{
	 U(r) &=& 4\epsilon\left(\left(\frac{\sigma}{r-\Delta}\right)^{12} - \left(\frac{\sigma}{r-\Delta}\right)^{6}\right) + U_{shift} & r < r_{cut} \\
	 &=& 0 & r \ge r_{cut}
	 \f}
	 */
	class slj : public Potential {
	public:
		void set_coeff (const vector <double> coeffs);
		inline double energy (const double r2);
		inline void force (Atom *a1, Atom *a2, const double r2, const double *xyz);	
		
	private:
		double epsilon_;
		double sigma_;
		double delta_;
	};	

	/*!
	 \param [in] \*coeffs Should be passed as [\epsilon, \sigma, \Delta]
	 */
	void slj::set_coeff (const vector <double> coeffs) {
		if (coeffs.size() != 3) {
			flag_error ("Failed to set coefficients for SLJ", __FILE__, __LINE__);
			exit(-1);
		}
		epsilon_ = coeffs[0];
		sigma_ = coeffs[1];
		delta_ = coeffs[2];
	}
	
	/*!
	 \param [in] \*a1 Pointer to atom1
	 \param [in] \*a2 Pointer to atom2
	 \param [in] r2 Minimum image distance squared between atoms
	 */
	inline double slj::energy (const double r2) {
		if (r2 < rcut*rcut) {
			double a = sigma_/(sqrt(r2)-delta_), a2, a6;
			a2 = a*a;
			a6 = a2*a2*a2;
			return 4.0*epsilon_*(a6*a6-a6)+Ushift;
		} else {
			return 0.0;
		}
	}
	
	/*!
	 The vector xyz MUST be pointing from a1 to a2 for the vectors to be correct.  This is automatically handled in min_image_dist2().
	 \param [in] \*a1 Pointer to atom1
	 \param [in] \*a2 Pointer to atom2
	 \param [in] r2 Minimum image distance squared between atoms
	 \param [in] \*xyz The minimum image vector from a1 to a2
	 */
	inline void slj::force (Atom *a1, Atom *a2, const double r2, const double *xyz) {
		if (r2 < rcut*rcut) {
			double r = sqrt(r2), b = 1.0/(r-delta_), a = sigma_*b, a2, a6, val, factor;
			a2 = a*a;
			a6 = a2*a2*a2;
			
			factor = 24.0*epsilon_*a6*(2.0*a6-1.0)*b/r;
			for (int i = 0; i < 3; ++i) {
				val = xyz[i]*factor;
				a1->force[i] -= val;
				a2->force[i] += val;
			}
		} else {
			return;
		}
	}
}

#endif