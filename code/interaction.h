/**
 Interaction Information
 \author Nathan A. Mahynski
 **/

#ifndef INTERACTION_H_
#define INTERACTION_H_

#include "atom.h"
#include "misc.h"

using namespace std;
using namespace atom;
using namespace misc;

// do these typedefs
typedef vector <double> (*force_ptr) (const double r2, const double *xyz, const vector <double> *args);
typedef double (*energy_ptr) (const double r2, const double *xyz, const vector <double> *args);

// declare new slj etc. functions here

void force_serial (Atom *a1, Atom *a2, const vector <double> *box);				//!< Does what inline vector <double> force () does, and apply the force to a1 and a2
	

//! This class stores how a pair of particles interacts
/*!
 Returns energy and force as long as r < r_{cut}, else 0.
 */
class Interaction {
public:
	vector <double> force (const Atom *a1, const Atom *a2, const vector <double> *box);	//!< Computes the cartesian force vector a1 experiences because of a2
	double energy (const Atom *a1, const Atom *a2, const vector <double> *box);			//!< Computes the energy between a1 and a2
	void set_force (force_ptr iforce) {my_force_ = iforce;}										//!< Assign the force calculator
	void set_energy (energy_ptr iener) {my_energy_ = iener;}									//!< Assign the potential calculator
	void set_args (const vector <double> args) {energy_args_ = args;}
	void set_force_args (const vector <double> args) {force_args_ = args;}
	void set_rcut (const double rcut) {rcut2_ = rcut_*rcut_;}
						  
private:
	double rcut2_;							//!< Note that rcut is stored as a square so that checking will be faster
	force_ptr my_force_;
	energy_ptr my_energy_;
	vector <double> force_args_;			//!< Force arguments
	vector <double> energy_args_;			//!< Energy arguments
};

//! Cosntructor sets U_{shift} and r_{cut} to 0.0 to will not compute anything until user sets these values.
Interaction::Interaction () {
	Ushift_ = 0.0;
	rcut_ = 0.0;
}


#endif


