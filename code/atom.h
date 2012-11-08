/**
 MD Atom Information
 \author Nathan A. Mahynski
**/

#ifndef ATOM_H_
#define ATOM_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <map>

using namespace std;

//! Namespace containing pertinent Atom information
namespace atom {
  
	class Atom {
	public:
		vector <double> pos() const {return pos_;}		//!< Report atom position
		vector <double> vel() const {return vel_;}		//!< Report atom velocity
		int type() const {return type_;}				//!< Report internally indexed atom type
		double mass() const {return mass_;}				//!< Report atomic mass
		double diam() const {return diam_;}				//!< Report atomic diameter
		int sys_index () const {return sys_index_;}		//!< Report index in System vector where this atom has been stored
		
		void set_pos (const vector<double> new_pos);			//!< Set the atom position 
		void set_vel (const vector<double> new_vel);			//!< Set the atom velocity
		void set_diam (const double diam) {diam_ = diam;}		//!< Set the atomic diameter
		void set_mass (const double mass) {mass_ = mass;}		//!< Set the atomic mass
		void set_type (const int type) {type_ = type;}			//!< Set the atomic type index
		void set_index (const int index) {sys_index_ = index;}	//!< Set the system index, which is the position in System's internal vector this atom is located
		
	private:
		vector <double> pos_;		//!< Cartesian coordinates
		vector <double> vel_;		//!< Cartesian velocities
		int type_;					//!< Internally indexed type of this atom
		double mass_;				//!< Atomic mass (in reduced units)
		double diam_;				//!< Atomic diameter (in reduced units)
		int	sys_index_;				//!< Location in internal vector of System::atoms_
	};
	
	// Can define addition "atom-relevant" functions, classes, etc. below...
	
}

#endif