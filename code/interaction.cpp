/**
 Interaction Information
 \authors{Nathan A. Mahynski, George A. Khoury}
 **/

#include "interaction.h"

class NewException : public exception {
	virtual const char* what() const throw() {
		return "Fene distance out of bounds (r > r0)";
	}
} fene_bounds_error;

// define new slj etc. functions here
/*void force_serial (Atom *atom1, Atom *atom2, const vector <double> *box) {
	double delta_=1., rcut=1., sigma_=0.1, epsilon_=0.1;
	// compute min image distance
	double xyz[3];
	double d2 = min_image_dist2 ((const Atom *)atom1, (const Atom *) atom2, box, xyz);
	
	double r = sqrt(d2), x = r - delta_;
	vector <double> force_vec(3,0.0);
	if (x < rcut) {
		double b = 1.0/x, a = sigma_*b, a2 = a*a, a6 = a2*a2*a2, val, factor;
		factor = 24.0*epsilon_*a6*(2.0*a6-1.0)*b/r;
		for (int i = 0; i < 3; ++i) {
			val = xyz[i]*factor;
			atom1->force[i] -= val;
			atom2->force[i] += val;
		}
	}
}*/

//! Shifted Lennard-Jones Force
/*!
 This is the same as standard LJ if \Delta = 0.  This is generally useful for systems with large size asymmetries. The energy U(r) is returned:
 \f{eqnarray*}{
 U(r) &=& 4\epsilon\left(\left(\frac{\sigma}{r-\Delta}\right)^{12} - \left(\frac{\sigma}{r-\Delta}\right)^{6}\right) + U_{shift}
 \f}
 Forces are stored on atoms: 
 \f[
 F_i = -\frac{\del U}{\del r}\frac{\del r}{\del x_i} = -\frac{\del U}{\del r}\frac{x_i}{r}
 \f]
 \param [in,out] \*atom1 Pointer to first atom
 \param [in,out] \*atom2 Pointer to second atom
 \param [in] \*box Pointer to vector of box size
 \param [in] \*args Vector of arguments <epsilon, sigma, delta, U_{shift}, rcut^2>
 */
double slj (Atom *atom1, Atom *atom2, const vector <double> *box, const vector <double> *args) {
	// compute min image distance
	double xyz[3];
	double d2 = min_image_dist2 ((const Atom *)atom1, (const Atom *) atom2, box, xyz);
	if (d2 < args->at(4)) {
		double delta=args->at(2), sigma=args->at(1), epsilon=args->at(0);
		double r = sqrt(d2), x = r - delta;	
		double b = 1.0/x, a = sigma*b, a2 = a*a, a6 = a2*a2*a2, val, factor;
		factor = 24.0*epsilon*a6*(2.0*a6-1.0)*b/r;
		// The vector xyz MUST be pointing from a1 to a2 for the vectors to be correct.  This is automatically handled in min_image_dist2().
		for (int i = 0; i < 3; ++i) {
			val = xyz[i]*factor;
			atom1->force[i] -= val;
			atom2->force[i] += val;
		}
		return 4.0*epsilon*(a6*a6-a6)+args->at(3);
	} else {
		return 0.0;
	}
}

//! Harmonic Bond
/*!
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
	double xyz[3];
	double d2 = min_image_dist2 ((const Atom *) a1, (const Atom *) a2, box, xyz);
	double d1 = sqrt(d2), factor = args->at(0)*(args->at(1)/d1-1.0);
	for (int i = 0; i < 3; ++i) {
		a1->force[i] -= xyz[i]*factor;
		a2->force[i] += xyz[i]*factor;
	}
	return 0.5*args->at(1)*(d1-args->at(0))*(d1-args->at(0));
}

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
 \param [in,out] \*atom1 Pointer to first atom
 \param [in,out] \*atom2 Pointer to second atom
 \param [in] \*box Pointer to vector of box size
 \param [in] \*args Vector of arguments <epsilon, sigma, delta, k, r0>
 */
double fene (Atom *a1, Atom *a2, const vector <double> *box, const vector <double> *args) {
	char err_msg[MYERR_FLAG_SIZE];
	double xyz[3];
	double d2 = min_image_dist2 ((const Atom *) a1, (const Atom *) a2, box, xyz);
	double d1 = sqrt(d2), d1shift = d1 - args->at(2);
	double factor = d1shift/(d1shift/args->at(4)*d1shift/args->at(4)-1.0)/d1, energy = 0.0;

        // check_fene
        // check if atom coordinates in a fene bond are further than r0 apart
        // this can cause a singularity so we want to guard against this
        if ( d1 > args->at(4) ) {
			throw(fene_bounds_error);
			return 0.0;
        }
	
	// compute logarithmic portion
	for (int i = 0; i < 3; ++i) {
		a1->force[i] -= xyz[i]*factor;
		a2->force[i] += xyz[i]*factor;
	}
	energy = -0.5*args->at(3)*args->at(4)*args->at(4)*log(1.0-d1shift/args->at(4)*d1shift/args->at(4));
	
	if (d1shift < WCA_CUTOFF*args->at(1)) {
		// use WCA portion of the potential
		double factor1 = (args->at(1)/d1shift), d2_wca = factor1*factor1, d6 = d2_wca*d2_wca*d2_wca;
		double factor2 = 24.0/d1shift*args->at(0)*d6*(2.0*d6-1.0)/d1;
		//		cout<<"d1 = "<<d1<<" , f = "<<factor<<" , f2 = "<<factor2<<endl;
		for (int i = 0; i < 3; ++i) {
			a1->force[i] -= xyz[i]*factor2;
			a2->force[i] += xyz[i]*factor2;
		}
		energy += (4.0*args->at(0)*d6*(d6-1.0)+args->at(0));
	} 
	
	return energy;
}
