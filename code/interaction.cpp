/**
 Interaction Information
 \author Nathan A. Mahynski
 **/

#include "interaction.h"

// define new slj etc. functions here
inline void Interaction::force_serial (Atom *atom1, Atom *atom2, const vector <double> *box) {
	double delta_=1., rcut=1., sigma_=0.1, epsilon_=0.1;
	// compute min image distance
	double xyz[3];
	double d2 = min_image_dist2 (atom1, atom2, box, xyz);
	
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
}

/*!
 \param [in] \*a1 Pointer to atom 1
 \param [in] \*a2 Pointer to atom 2
 \param [in] \*box Pointer to vector of cartesian box size
*/
inline vector <double> Interaction::force (const Atom *a1, const Atom *a2, const vector <double> *box) {
	// compute min image distance
	double xyz[3];
	double d2 = min_image_dist2 (a1, a2, box, xyz);
	if (d2 < rcut2_) {
		return my_force_ (d2, xyz, &force_args_);
	} else {
		vector <double> noforce(3,0.0);
		return noforce;
	}
}

/*!
 \param [in] \*a1 Pointer to atom 1
 \param [in] \*a2 Pointer to atom 2
 \param [in] \*box Pointer to vector of cartesian box size
 */
inline double Interaction::energy (const Atom *a1, const Atom *a2, const vector <double> *box) {
	// compute min image distance
	double xyz[3];
	double d2 = min_image_dist2 (a1, a2, box, xyz);
	if (d2 < rcut2_) {
		return my_energy_ (d2, xyz, &energy_args_);
	} else {
		return 0.0;
	}
}

//! Shifted Lennard-Jones Force
/*!
 This is the same as standard LJ if \Delta = 0.  This is generally useful for systems with large size asymmetries.
 \f{eqnarray*}{
 U(r) &=& 4\epsilon\left(\left(\frac{\sigma}{r-\Delta}\right)^{12} - \left(\frac{\sigma}{r-\Delta}\right)^{6}\right) + U_{shift}
 \f}
 Returns 
 \f[
 F_i = -\frac{\del U}{\del r}\frac{\del r}{\del x_i} = -\frac{\del U}{\del r}\frac{x_i}{r}
 \f]
 The vector xyz MUST be pointing from a1 to a2 for the vectors to be correct.  This is automatically handled in min_image_dist2().
 \param [in] r2 Minimum image distance squared between atoms
 \param [in] \*xyz The minimum image vector from a1 to a2
 \param [in] \*args Vector of arguments <epsilon, sigma, delta>
 */
vector<double> slj_force (const double r2, const double *xyz, const vector <double> *args) {
	double r = sqrt(r2), x = r - args->at(2);
	vector <double> force_vec(3,0.0);
	
	// checking if r < rcut is done beforehand so do not need it here
	double b = 1.0/x, a = args->at(1)*b, a2 = a*a, a6 = a2*a2*a2, val, factor;
	factor = 24.0*args->at(0)*a6*(2.0*a6-1.0)*b/r;
	for (int i = 0; i < 3; ++i) {
		force_vec[i] = xyz[i]*factor;
	}
	return force_vec;
}

//! Shifted Lennard-Jones Energy
/*!
 This is the same as standard LJ if \Delta = 0.  This is generally useful for systems with large size asymmetries.
 \f{eqnarray*}{
 U(r) &=& 4\epsilon\left(\left(\frac{\sigma}{r-\Delta}\right)^{12} - \left(\frac{\sigma}{r-\Delta}\right)^{6}\right) + U_{shift} 
 \f}
 \param [in] r2 Minimum image distance squared between atoms
 \param [in] \*xyz The minimum image vector from a1 to a2
 \param [in] \*args Vector of arguments <epsilon, sigma, delta, U_{shift}>
 */
double slj_energy (const double r2, const double *xyz, const vector <double> *args) {
	double r = sqrt(r2), x = r - args->at(2);
	
	// checking if r < rcut is done beforehand so do not need it here
	double a = args->at(1)/x, a2 = a*a, a6 = a2*a2*a2;
	return 4.0*args->at(0)*(a6*a6-a6)+args->at(3);
}


