/**
 Interaction Information
 \author Nathan A. Mahynski
 **/

#ifndef INTERACTION_H_
#define INTERACTION_H_

#include "atom.h"
#include "misc.h"
#include "global.h"

using namespace std;
using namespace atom;
using namespace misc;

// Functor for functions that compute (and store) force vector between 2 atoms, and return the energy between them.
typedef double (*force_energy_ptr) (Atom *a1, Atom *a2, const vector <double> *box, const vector <double> *args);

double slj (Atom *a1, Atom *a2, const vector <double> *box, const vector <double> *args);						//!< Computes force and energy of Shifted Lennard-Jones interaction
double fene (Atom *a1, Atom *a2, const vector <double> *box, const vector <double> *args);						//!< Computes force and energy of a Fene bond
double harmonic (Atom *a1, Atom *a2, const vector <double> *box, const vector <double> *args);					//!< Computes force and energy of a Harmonic bond
	
//! This class stores how a pair of particles interacts
/*!
 Returns energy and force as long as r < r_{cut}, else 0.
 */
class Interaction { 
public:
	Interaction() {my_force_energy_ = NULL;} 
	~Interaction() {};
	double force_energy (Atom *a1, Atom *a2, const vector <double> *box) {return my_force_energy_ (a1, a2, box, &energy_args_);}					//!< Computes force (stored on atoms) and energy (returned)
	void set_force_energy (force_energy_ptr ife) {my_force_energy_ = ife;}									//!< Assign the potential calculator
	void set_args (const vector <double> args) {energy_args_ = args;}
	force_energy_ptr check_force_energy_function () const {return my_force_energy_;}						//!< Return the function for force and energy calculations
						  
private:
	force_energy_ptr my_force_energy_;
	vector <double> energy_args_;			//!< Energy (and force) arguments
};

#endif


