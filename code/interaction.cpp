/**
 Interaction Information
 \author Nathan A. Mahynski
 **/

#include "interaction.h"

// define new slj etc. functions here
void force_serial (Atom *atom1, Atom *atom2, const vector <double> *box) {
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
	return my_force_ (d2, xyz, &force_args_);
}

/*!
 \param [in] \*a1 Pointer to atom 1
 \param [in] \*a2 Pointer to atom 2
 \param [in] \*box Pointer to vector of cartesian box size
 */
inline double Interaction::energy(const Atom *a1, const Atom *a2, const vector <double> *box) {
	// compute min image distance
	double xyz[3];
	double d2 = min_image_dist2 (a1, a2, box, xyz);
	return my_energy_ (d2, xyz, &energy_args_);
}



