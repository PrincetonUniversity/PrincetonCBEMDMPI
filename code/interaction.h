/**
 Interaction Information
 \author Nathan A. Mahynski
 **/

#ifndef INTERACTION_H_
#define INTERACTION_H_

#include "bond.h"
#include "pair_potential.h"

using namespace std;
using namespace misc;
using namespace bond;
using namespace pair_potential;


// do these typedefs
typedef force_ptr vector <double> (*force_function) (const double r2, const double *xyz, const vector <double> *args);
typedef energy_ptr double (*energy_function) (const double r2, const double *xyz, const vector <double> *args);


// define new slj etc. functions here








//! This class stores how a pair of particles interacts
class Interaction {
public:
	inline vector <double> force (const Atom *a1, const Atom *a2, const vector <double> *box);	//!< Computes the cartesian force vector a1 experiences because of a2
	inline double energy (const Atom *a1, const Atom *a2, const vector <double> *box);			//!< Computes the energy between a1 and a2
	void set_force (force_ptr iforce) {my_force_ = iforce;}										//!< Assign the force calculator
	void set_energy (energy_ptr iener) {my_energy_ = iener;}									//!< Assign the potential calculator
	
private:
	force_ptr my_force_;
	energy_ptr my_energy_;
	vector <double> force_args_;			//!< Force arguments
	vector <double> energy_args_;			//!< Energy arguments
};

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


#endif