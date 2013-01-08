/*!
 \brief Source code for interaction functions
 \file interaction.cpp
 \authors{Nathan A. Mahynski, George Khoury}
 **/

#include "interaction.h"

//! Shifted Lennard-Jones Force
/*!
 This is the same as standard LJ if \Delta = 0.  This is generally useful for systems with large size asymmetries. The energy U(r) is returned:
 \f[
 U(r) &=& 4\epsilon\left(\left(\frac{\sigma}{r-\Delta}\right)^{12} - \left(\frac{\sigma}{r-\Delta}\right)^{6}\right) + U_{shift} & r - \Delta < r_{cut}  &=& 0 & r - \Delta \ge r_{cut}
 \f]
 Forces are added to atoms: 
 \f[
 F_i &=& -\frac{\del U}{\del r}\frac{\del r}{\del x_i} = -\frac{\del U}{\del r}\frac{x_i}{r} & r - \Delta < r_{cut} &=& 0 & r - \Delta \ge r_{cut}
 \f]
 \param [in,out] \*atom1 Pointer to first atom
 \param [in,out] \*atom2 Pointer to second atom
 \param [in] \*box Pointer to vector of box size
 \param [in] \*args Pointer to vector of arguments <epsilon, sigma, delta, U_{shift}, rcut^2>
 */
double slj (Atom *atom1, Atom *atom2, const vector <double> *box, const vector <double> *args) {
	// Compute min image distance
	double xyz[NDIM];
	double d2 = min_image_dist2 ((const Atom *)atom1, (const Atom *) atom2, box, xyz);
	double delta=args->at(2), r = sqrt(d2), x = r - delta;	
	
	// Check that r > delta, else error has occurred
	if (x < 0) {
		SljException slj_bounds_error (atom1->sys_index, atom2->sys_index, d2, args->at(2));
		throw(slj_bounds_error);
		return 0.0;
	}
	
	// If (r-delta)^2 < rcut^2 compute
	if (x*x < args->at(4)) {
		double sigma=args->at(1), epsilon=args->at(0);
		double b = 1.0/x, a = sigma*b, a2 = a*a, a6 = a2*a2*a2, val, factor;
		factor = 24.0*epsilon*a6*(2.0*a6-1.0)*b/r;
		/* The vector xyz MUST be pointing from a1 to a2 for the vectors to 
        be correct.  This is automatically handled in min_image_dist2().
        */
		for (int i = 0; i < NDIM; ++i) {
			val = xyz[i]*factor;
			atom1->force[i] -= val;
			atom2->force[i] += val;
		}
		return 4.0*epsilon*(a6*a6-a6)+args->at(3);
	} else {
		return 0.0;
	}
}

/*!
 Harmonic Bond
 The Harmonic bond potential is given by:
 \f[
 U(r) = \frac{1}{2}k\left(r - r_{0}\right)^2
 \f]
 \param [in,out] \*atom1 Pointer to first atom
 \param [in,out] \*atom2 Pointer to second atom
 \param [in] \*box Pointer to vector of box size
 \param [in] \*args Vector of arguments <k, r0>
 */
double harmonic (Atom *a1, Atom *a2, const vector <double> *box, const vector <double> *args) {
	double xyz[NDIM];
	double d2 = min_image_dist2 ((const Atom *) a1, (const Atom *) a2, box, xyz);
	double d1 = sqrt(d2), factor = args->at(0)*(1.0-args->at(1)/d1);
	for (int i = 0; i < NDIM; ++i) {
		a1->force[i] -= xyz[i]*factor;
		a2->force[i] += xyz[i]*factor;
	}
	return 0.5*args->at(0)*(d1-args->at(1))*(d1-args->at(1));
}

/*!
 Finitely Extensible Non-linear Elastic Bond (FENE)
 The Fene bond potential is given by:
 \f[
 U(r) = -\frac{1}{2}kr_{0}^2\text{ln}\left(1-\left(\frac{r-\Delta}{r_0}\right)^2\right) + U_{WCA}
 \f]
 Where the short range repulsion is provided by the WCA potential:
 \f[
 U_{WCA} &=& 4\epsilon \left( \left( \frac{\sigma}{r-\Delta} \right)^{12} - \left( \frac{\sigma}{r-\Delta} \right)^6 \right) + \epsilon & r < 2^{1/6}\sigma+\Delta &=& 0 &r-\Delta \ge 2^{1/6}\sigma
 \f]
 \Delta is usually set such that $\Delta = (d_i+d_j)/2-1$ where $d_i$ is the diameter of species i, but the user may decide on other parameters.
 \param [in,out] \*atom1 Pointer to first atom
 \param [in,out] \*atom2 Pointer to second atom
 \param [in] \*box Pointer to vector of box size
 \param [in] \*args Vector of arguments <epsilon, sigma, delta, k, r0>
 */
double fene (Atom *a1, Atom *a2, const vector <double> *box, const vector <double> *args) {
	double xyz[NDIM];
	double d2 = min_image_dist2 ((const Atom *) a1, (const Atom *) a2, box, xyz);
	double d1 = sqrt(d2), d1shift = d1 - args->at(2);
	double factor = args->at(3)*d1shift/(d1shift/args->at(4)*d1shift/args->at(4)-1.0)/d1, energy = 0.0;

	// Check if atom coordinates in a fene bond are further than r0 apart; this can cause a singularity so we want to guard against this.
	if (d1 > args->at(4)) {
		FeneException fene_bounds_error (a1->sys_index, a2->sys_index, d1, args->at(4));
		throw(fene_bounds_error);
		return 0.0;
	}
	
	// WCA portion potentially also has singularity that needs to be checked
	if (d1shift < 0) {
		SljException slj_bounds_error (a1->sys_index, a2->sys_index, d1, args->at(2));
		throw(slj_bounds_error);
		return 0.0;
	}
	
	// Compute logarithmic portion
	for (int i = 0; i < NDIM; ++i) {
		a1->force[i] -= xyz[i]*factor;
		a2->force[i] += xyz[i]*factor;
	}
	energy = -0.5*args->at(3)*args->at(4)*args->at(4)*log(1.0-d1shift/args->at(4)*d1shift/args->at(4));
	
	if (d1shift < WCA_CUTOFF*args->at(1)) {
		// use WCA portion of the potential
		double factor1 = (args->at(1)/d1shift), d2_wca = factor1*factor1, d6 = d2_wca*d2_wca*d2_wca;
		double factor2 = 24.0/d1shift*args->at(0)*d6*(2.0*d6-1.0)/d1;
		for (int i = 0; i < NDIM; ++i) {
			a1->force[i] -= xyz[i]*factor2;
			a2->force[i] += xyz[i]*factor2;
		}
		energy += (4.0*args->at(0)*d6*(d6-1.0)+args->at(0));
	} 
	
	return energy;
}
