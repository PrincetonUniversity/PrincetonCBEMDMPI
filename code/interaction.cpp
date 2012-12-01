/**
 Interaction Information
 \author Nathan A. Mahynski
 **/



#include "interaction.h"

// define new slj etc. functions here


/*!
 \param [in] \*a1 Pointer to atom 1
 \param [in] \*a2 Pointer to atom 2
 \param [in] \*box Pointer to vector of cartesian box size
*/
inline vector <double> Interaction::force (const Atom *a1, const Atom *a2, const vector <double> *box) {
	// compute min image distance
	double xyz[3];
	double d2 = min_image_dist2 (a1, a2, box, xyz);
	return my_force_ (d2, xyz, force_args_);
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
	return my_energy_ (d2, xyz, energy_args_);
}



