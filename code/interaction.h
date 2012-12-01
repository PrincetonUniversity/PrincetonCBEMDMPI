/*
 *  interaction.h
 *  
 *
 *  Created by Nathan A. Mahynski on 11/30/12.
 *  Copyright 2012 Princeton University. All rights reserved.
 *
 */

#ifndef INTERACTION_H_
#define INTERACTION_H_

#include "bond.h"
#include "pair_potential.h"

using namespace std;
using namespace misc;
using namespace bond;
using namespace pair_potential;

//! This class stores how a pair of particles interacts
class Interaction {
public:
	inline vector <double> force (const Atom *a1, const Atom *a2, const vector <double> *box);	//!< Computes the cartesian force vector a1 experiences because of a2
	inline double energy (const Atom *a1, const Atom *a2, const vector <double> *box);			//!< Computes the energy between a1 and a2
	void set_bond (Bond ibond) {bond_interaction_ = ibond; is_bonded_ = 1;}						//!< Assign the bonding potential
	void set_ppot (Potential ipair) {ppot_ = ipair; is_bonded_ = 0;}							//!< Assign the pair potential
	
private:
	int is_bonded_;			//!< If this is 1, the atoms are bonded so use bond potential, else use pair potential
	Bond bond_interaction_;	//!< Bond class
	Potential ppot_;		//!< Pair potential class
};

/*!
 \param [in] \*a1 Pointer to atom 1
 \param [in] \*a2 Pointer to atom 2
 \param [in] \*box Pointer to vector of cartesian box size
*/
inline vector <double> Interaction::force(const Atom *a1, const Atom *a2, const vector <double> *box) {
	// compute min image distance
	double xyz[3];
	double d2 = min_image_dist2 (a1, a2, box, xyz);
	
	// rcut belongs to these classes and they will establish if the force needs to be computed
	if (is_bonded_) {
		return bond_interaction_.force(d2, xyz);
	} else {
		return ppot_.force(d2, xyz);
	}
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
	
	// rcut belongs to these classes and they will establish if the force needs to be computed
	if (is_bonded_) {
		return bond_interaction_.energy(d2);
	} else {
		return ppot_.energy(d2);
	}
}


#endif