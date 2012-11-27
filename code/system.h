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
#include <list>
#include <algorithm>

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
		int set_box (const double* new_box);			//!< Set the global system box size
		void set_T (const double T) {Temp_ = T;}		//!< Set the system temperature
		void set_P (const double P) {Press_ = P;}		//!< Set the system pressure
		int natoms () const {return atoms_.size();}		//!< Return the number of atoms currently in this system (processor)
		int add_atom_type (const string atom_name);		//!< Index an atom name
		int add_bond_type (const string bond_name);		//!< Index a bond name
		int add_ppot_type (const string ppot_name);		//!< Index a pair potential's name
		double T() const {return Temp_;}				//!< Report the temperature of the system
		double P() const {return Press_;}				//!< Report the pressure of the system
		vector <double> box() const {return box_;}		//!< Report system size
		int add_atoms (Atom *new_atoms);				//!< Add atom(s) to the system 
		int delete_atoms (const vector <int> indices);	//!< Pop atoms with local indices from local storage
		inline Potential* ppot (const int type1, const int type2) const {return pmatrix_[type1][type2];}
		inline Potential* ppot (const Atom* a1, const Atom* a2) const {return pmatrix_[a1->type][a2->type];}
			
	private:
		double Temp_;									//!< System temperature in reduced units (kT)
		double Press_;									//!< System pressure in reduced units
		double KE_, U_;									//!< Total internal kinetic energy and potential energy
		vector <double> box_;							//!< Global system cartesian dimensions
		vector <Atom> atoms_;							//!< Vector of Atoms in the system
		vector <Bond> bond_;							//!< Vector of bond types in the system
		vector < vector < *Potential> > pmatrix_;		//!< Matrix of pair potential pointers based on internal index
		vector <double> masses_;						//!< Vector containing masses of each type of Atom in the system
		unordered_map <string, int> atom_type_;			//!< Maps user specified name of atom type to internal index
		unordered_map <string, int> bond_type_;			//!< Maps user specified name of bond type to internal index
		unordered_map <string, int> ppot_type_;			//!< Maps user specified name of pair potential type to an internal index
		unordered_map <int, int> glob_to_loc_id_;		//!< Maps global sys_index to the local index of atoms_ an atom is stored at on each processor; the opposite conversion can be done with lookup of Atom::sys_index
		vector <int> on_proc_;							//!< Stores processor each atom (by global sys_index) "lives" on
		int num_neighbor_procs;							//!< Number of neighboring processors this system must communicate with
		int *neighbor_procs_;							//!< Array of ranks of neighboring processors this system must communicate with
	};
	
	//! Initializes a system from a file, automatically handles MPI
	int initialize (const char *filename, int *rank, int *nprocs, System *sys);
	
	//!< Finalizes a system by performing any last minute tasks that must be done before the simulation finishes.
	int finalize();
	
	//!< Decomposes a box into domains for each processor to handle
	int domain_decomposition (const vector <double> box, const int nprocs, vector < vector <int> > *neighbors);
	
	//! Read atom coordinates and properties from a file (xml)
	int read_xml (const string filename, const int nprocs, System *new_systems);			
	
	// declare system before MPI Init
	// then after initialization each proc has a copy of an empty System
	// then read from file from proc1 and SEND data to each as it is read in!
	// i.e. things like T/P are set for all, but atoms only need to be sent to the necessary dest but with global id's
	// collect and create atoms before hand so that bonds can be assigned before being sent out
	
	// bonding -- store global index of atoms each is bonded to on the Atom class so bonds for atoms on each each proc can be cycled quickly
	
	// for now do simple force compute on each proc
	// will need to establish r_cut,max for this
	// at the beginning of this routine check for atoms close to boundary and send an array as necessary to neighbors
	// then loop atoms that belong ON that proc with others on same proc and those received to compute their forces
	
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
	 Remove atoms from the system.  Needs to sort indices because erase() operation reorders things; also, because of this it is fastest to pop from lowest to highest index.
	 Returns the number of atoms deleted.
	 */
	int delete_atoms (const vector <int> indices) {
		vector <Atom>::iterator it = atoms_.begin();
		
		// sort indices from lowest to highest
		sort (indices.begin(), indices.end());
		
		// pop in this order
		int shift = 0;
		for (int i = 0; i < indices.size(); ++i) {
			// update glob_to_loc for all atoms following the erased atoms
			if (i < indices.size()-1) {
				upper = indices[i+1]-shift;
			} else {
				upper = atoms_.size();
			}
			for (int j = indices[i]+1-shift; j < upper; ++j) {
				glob_to_loc[atoms_[j].sys_index] -= (shift+1);
			}

			atoms_.erase(it-shift+indices[i]);
			++shift;
		}
		return shift;
	}
				 
	/*!
	 Attempt to push an atom(s) into the system.  This assigns the map automatically to link the atoms global index to the local storage location.
	 This reallocates the internal vector that stores the atoms; if a memory error occurs during such reallocation, an error is given and the system exits.
	 \param [in] natoms Length of the array of atoms to add to the system.
	 \param [in] \*new_atoms Pointer to an array of atoms the user has created elsewhere.
	 \param [out] \*update_proc Array of global indices that have just been added to this proc to be used for updating on_proc_ list.
	 */
	int* System::add_atoms (const int natoms, Atom *new_atoms) {
		int index = atoms_.size(), *update_proc = new int [natoms];
		if (natoms < 1) {
			return natoms;
		}
		
		for (int i = 0; i < natoms; ++i) {
			try {
				atoms_.push_back(new_atoms[i]);
			}
			catch (bad_alloc& ba) {
				sprintf(err_msg, "Could not allocate space for new atoms in the system");
				flag_error (err_msg, __FILE__, __LINE__);
				finalize();
				exit(BAD_MEM);
			}
			glob_to_loc_id_[new_atoms[i].sys_index] = index;
			update_proc[i] = new_atoms[i].sys_index;
			index++;
		}

		return update_proc;
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
	
	// Can define addition "system-relevant" functions, classes, etc. below...

}

#endif