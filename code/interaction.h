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

//! Fene exception class is thrown if there is an error
class FeneException : public exception {
public:
	FeneException (const int ind1, const int ind2, const double dist, const double r0) {ind1_=ind1; ind2_=ind2; dist_=dist, r0_=r0;}
	virtual const char* what() const throw() {
		char err_msg[1000];
		sprintf(err_msg, "Fene distance for atoms (%d,%d) = %g > %g, out of bounds (r > r0)", ind1_, ind2_, dist_, r0_);
		return err_msg;
	}
protected:
	int ind1_, ind2_;
	double dist_, r0_;
};

//! SLJ exception class is thrown if there is an error
class SljException : public exception {
public:
	SljException (const int ind1, const int ind2, const double dist, const double delta) {ind1_=ind1; ind2_=ind2; dist_=dist, delta_=delta;}
	virtual const char* what() const throw() {
		char err_msg[1000];
		sprintf(err_msg, "Shifted Lennard-Jones for atoms (%d,%d) has separation of %g < %g out of bounds (r < delta)", ind1_, ind2_, dist_, delta_);
		return err_msg;
	}
protected:
	int ind1_, ind2_;
	double dist_, delta_;
};

#endif


