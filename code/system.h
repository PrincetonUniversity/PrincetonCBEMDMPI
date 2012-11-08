/**
 MD System Information
 \author Nathan A. Mahynski
**/

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include "atom.h"
#include "bond.h"
#include "misc.h"

using namespace std;
using namespace atom;
using namespace bond;
using namespace misc;

//! Namespace containing all system variable and quantities
namespace system {
	char err_msg[ERR_FLAG_SIZE]; //!< Error message buffer commonly used in routines in this namespace
	
	class System {
	public:
		System();
		~System();
		int read_xml (const string filename);			//!< Read atom coordinates and properties from a file (xml)
		int set_integrator (const integrator new_int);	//!< Set the system integrator
		int set_box (const double* new_box);			//!< Set the system box size
		void set_T (const double T) {Temp_ = T;}		//!< Set the system temperature
		void set_P (const double P) {Press_ = P;}		//!< Set the system pressure
		int integrate ();								//!< Integrate the system forward in time one step
		int natoms () const {return atoms_.size();}		//!< Return the number of atoms currently in the system
		int nbonds() const {return bonds_.size();}		//!< Return the number of bonds currently in the system
		int add_atom_type (const string atom_name);		//!< Index an atom name
		int add_bond_type (const string bond_name);		//!< Index a bond name
		int add_ppot_type (const string ppot_name);		//!< Index a pair potential's name
		double T() const {return Temp_;}				//!< Report the temperature of the system
		double P() const {return Press_;}				//!< Report the pressure of the system
		vector <double> box() const {return box_;}		//!< Report system size
		double dt() const {return dt_;}					//!< Report the timestep being used by the system
		int add_atoms (const vector <Atom> *new_atoms);	//!< Add atoms to the system -- this should update the atoms as well with indices
		int add_bond (const int atom1, const int atom2, const string type);		//!< Add a bond between atoms

	private:
		double Temp_;									//!< System temperature in reduced units (kT)
		double Press_;									//!< System pressure in reduced units
		vector <double> box_;							//!< System cartesian dimensions
		double dt_;										//!< Timestep size
		vector <Atom> atoms_;							//!< Vector of Atoms in the system
		vector < pair < pair <int, int>, int> > bonded_;//!< Vector of indices of atoms in atoms_ that are bonded, and the type of bond they have
		vector <Bond> bond_;							//!< Vector of bond types in the system
		vector <double> masses_;						//!< Vector containing masses of each type of Atom in the system
		Integrator sys_integrator_;						//!< System integrator function
		map <string, int> atom_type_;					//!< Maps user specified name of atom type to internal index
		map <string, int> bond_type_;					//!< Maps user specified name of bond type to internal index
		map <string, int> ppot_type_;					//!< Maps user specified name of pair potential type to an internal index
	};
	
	/*!
	 Upon initialization, resize vectors as necessary.  
	 */
	System::System() {
		try {
			box_.resize(3,-1);
		}
		catch (bad_alloc& ba) {
			sprintf(err_msg, "Could not allocate space for system's box size vector");
			flag_error (err_msg, __FILE__, __LINE__);
			exit(BAD_MEM);
		}
	}
				 
	/*!
	 Attempt to push an atom(s) into the system, these are stored in a vector and the index at which it is stored can be reported using Atom::sys_index().
	 This function returns an integer value of the number of atoms added to the system so the user can check that number was what they wanted. This 
	 reallocates the internal vector that stores the atoms; if a memory error occurs during such reallocation, an error is given and the system exits.
	 \param [in] \*new_atoms Pointer to a vector of atoms the user has created elsewhere; a copy is pushed into the system leaving the originals intact.
	 */
	int System::add_atoms (const vector <Atom> *new_atoms) {
		int natoms = new_atoms->size(), index = atoms_.size();
		if (natoms < 1) {
			return natoms;
		}
		
		for (int i = 0; i < natoms; ++i) {
			try {
				atoms_.push_back(new_atoms->at(i));
			}
			catch (bad_alloc& ba) {
				sprintf(err_msg, "Could not allocate space for new atoms in the system");
				flag_error (err_msg, __FILE__, __LINE__);
				exit(BAD_MEM);
			}
			atoms_[index].set_index(index);
			index++;
		}
		
		return natoms;
	}
	
	/*!
	 Adds a bond between two atoms.
	 */
	int System::add_bond (const int atom1, const int atom2, const string type) {
		;
	}
	
	/*!
	 Tries to add an atom type to the system, associating a user specified name with an internal index to reference this type in the future.
	 This can return 3 different values:
	 \par
	 Returns 0 if, atom_name was new and was successfully indexed.
	 \par
	 Returns -1 if, atom_name was bad (i.e. empty string).
	 \par
	 Returns +1 if, atom_name already exists in the system and could not be indexed again.
	 \par
	 \param [in] atom_name User specified name to index
	 */
	int System::add_atom_type (const string atom_name) {
		if (atom_name.size() == 0) return -1;
		if (atom_type_.find(atom_name) == atom_type_.end()) {
			int size = atom_type_.size();
			atom_type_[atom_name] = size;
			return 0;
		} else {
			return 1;
		}
	}
	
	/*!
	 Tries to add a bond type to the system, associating a user specified name with an internal index to reference this type in the future.
	 This can return 3 different values:
	 \par
	 Returns 0 if, bond_name was new and was successfully indexed.
	 \par
	 Returns -1 if, bond_name was bad (i.e. empty string).
	 \par
	 Returns +1 if, bond_name already exists in the system and could not be indexed again.
	 \par
	 \param [in] bond_name User specified name to index
	 */
	int System::add_bond_type (const string bond_name) {
		if (bond_name.size() == 0) return -1;
		if (bond_name.find(bond_name) == bond_type_.end()) {
			int size = bond_name.size();
			bond_type_[bond_name] = size;
			return 0;
		} else {
			return 1;
		}
	}
	
	/*!
	 Tries to add a pair potential type to the system, associating a user specified name with an internal index to reference this type in the future.
	 This can return 3 different values:
	 \par
	 Returns 0 if, ppot_name was new and was successfully indexed.
	 \par
	 Returns -1 if, ppot_name was bad (i.e. empty string).
	 \par
	 Returns +1 if, ppot_name already exists in the system and could not be indexed again.
	 \par
	 \param [in] ppot_name User specified name to index
	 */
	int System::add_ppot_type (const string ppot_name) {
		if (ppot_name.size() == 0) return -1;
		if (ppot_name.find(ppot_name) == ppot_type_.end()) {
			int size = ppot_name.size();
			ppot_type_[ppot_name] = size;
			return 0;
		} else {
			return 1;
		}
	}

	//!< Run (i.e. integrate) a system forward in time for a specified number of timesteps
	/*!
	 Function returns 0 if successful, -1 if it encountered an error.
	 \param [in] timesteps Number of timesteps to integrate over
	 */
	int run (const int timeteps) {
		// The way this function is written it can be easily interpreted by SWIG with python!
		// Inside this function we should initially check if the system is "prepared", i.e. all necessary vars specified
		
		// loop for timestep steps, report any errors as encountered
		
		// That's it!
	}
	
	// Can define addition "system-relevant" functions, classes, etc. below...
	
	
	
	
	
	
	
	
	// need integrator class (make it inheritable so that sub-classes of integrators can be defined)

}

#endif