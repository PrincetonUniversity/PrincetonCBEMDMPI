/**
 Interaction Information
 \author Nathan A. Mahynski
 **/

#ifndef INTERACTION_H_
#define INTERACTION_H_

#include "atom.h"
#include "bond.h"
#include "misc.h"
#include "pair_potential.h"

using namespace std;
using namespace atom;
using namespace misc;
using namespace bond;
using namespace pair_potential;


// do these typedefs
typedef vector <double> (*force_ptr) (const double r2, const double *xyz, const vector <double> *args);
typedef double (*energy_ptr) (const double r2, const double *xyz, const vector <double> *args);


// declare new slj etc. functions here



//! This class stores how a pair of particles interacts
class Interaction {
public:
	inline vector <double> force (const Atom *a1, const Atom *a2, const vector <double> *box);	//!< Computes the cartesian force vector a1 experiences because of a2
	inline void force_serial (Atom *a1, Atom *a2, const vector <double> *box);				//!< Does what inline vector <double> force () does, and apply the force to a1 and a2
	inline double energy (const Atom *a1, const Atom *a2, const vector <double> *box);			//!< Computes the energy between a1 and a2
	void set_force (force_ptr iforce) {my_force_ = iforce;}										//!< Assign the force calculator
	void set_energy (energy_ptr iener) {my_energy_ = iener;}									//!< Assign the potential calculator
	
private:
	force_ptr my_force_;
	energy_ptr my_energy_;
	vector <double> force_args_;			//!< Force arguments
	vector <double> energy_args_;			//!< Energy arguments
};


#endif


