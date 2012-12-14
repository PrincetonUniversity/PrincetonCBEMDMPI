/**
MD System Information
\author Nathan A. Mahynski
**/

#include "system.h"

namespace sim_system {	
	/*!
	 Upon initialization, resize vectors as necessary.  
	 */
	System::System() {
		try {
			box_.resize(3,-1);
		}
		catch (bad_alloc& ba) {
			char err_msg[MYERR_FLAG_SIZE]; 
			sprintf(err_msg, "Could not allocate space for system's box size vector");
			flag_error (err_msg, __FILE__, __LINE__);
			exit(BAD_MEM);
		}
	}

	//! Default destructor
	System::~System() {
		;
	}


	/*!
	 Remove atoms from the system.  Needs to sort indices because erase() operation reorders things; also, because of this it is fastest to pop from lowest to highest index.
	 Returns the number of atoms deleted.
	 */
	int System::delete_atoms (vector <int> indices) {
		vector <Atom>::iterator it = atoms_.begin();
  
		// sort indices from lowest to highest
		sort (indices.begin(), indices.end());
  
		// pop in this order
		int shift = 0;
		int upper;
		for (int i = 0; i < indices.size(); ++i) {
			// update glob_to_loc for all atoms following the erased atoms
			if (i < indices.size()-1) {
				upper = indices[i+1]-shift;
			} else {
				upper = atoms_.size();
			}
			for (int j = indices[i]+1-shift; j < upper; ++j) {
				glob_to_loc_id_[atoms_[j].sys_index] -= (shift+1);
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
	 \param [out] \*update_proc Array of global indices that have just been added to this proc.
	 */
	int* System::add_atoms (const int natoms, Atom *new_atoms) {
		int index = atoms_.size(), *update_proc = new int [natoms];
		for (int i = 0; i < natoms; ++i) {
			try {
				atoms_.push_back(new_atoms[i]);
			}
			catch (bad_alloc& ba) {
				char err_msg[MYERR_FLAG_SIZE]; 
				sprintf(err_msg, "Could not allocate space for new atoms in the system");
				flag_error (err_msg, __FILE__, __LINE__);
				//finalize();
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
		if (bond_type_.find(bond_name) == bond_type_.end()) {
			int size = bond_name.size();
			bond_type_[bond_name] = size;
			return 0;
		} else {
			return 1;
		}
	}

	
	/*!
	 Returns the internal index associated with this atom name; returns -1 if not found.
	 \param [in] name User defined name of atom type.
	 */
	inline int System::atom_type (const string name) {
		if (name.size() == 0) return -1;
		if (atom_type_.find(name) == atom_type_.end()) {
			return -1;
		} else {
			return atom_type_.find(name)->second;
		}
	}
	
	/*!
	 Returns the internal index associated with this bond name; returns -1 if not found.
	 \param [in] name User defined name of bond type.
	 */
	inline int System::bond_type (const string name) {
		if (name.size() == 0) return -1;
		if (bond_type_.find(name) == bond_type_.end()) {
			return -1;
		} else {
			return bond_type_.find(name)->second;
		}
	}

	/*!
	 Adds a new bond and associated information to the System object.
	 \param [in] atom1 Global index of atom1 in the bond.
	 \param [in] atom1 Global index of atom2 in the bond.
	 \param [in] type Internal index associated with this bond type.
	 */
	inline void System::add_bond (const int atom1, const int atom2, const int type) {
		pair <int, int> new_bond = make_pair (atom1, atom2);
		bonded_.push_back(new_bond);
		bonded_type_.push_back(type);
	}
	
	/*!
	 Returns the string "NULL" if failed, else user defined name of atom.
	 \param [in] index Internal index to locate and return the name associated.
	 */
	inline string System::atom_name (const int index) {
		string name = "NULL";
		char err_msg[MYERR_FLAG_SIZE];
		if (index >= atom_type_.size()) {
			sprintf(err_msg, "Index %d out of range", index);
			flag_error (err_msg, __FILE__, __LINE__);
			return name;
		} else {
			typedef std::map<std::string, int>::iterator it_type;
			for	(it_type iterator = atom_type_.begin(); iterator != atom_type_.end(); iterator++) {
				if (iterator->second == index) {
					return iterator->first;
				}
			}
			sprintf(err_msg, "Could not locate index %d in indexed atom types", index);
			flag_error (err_msg, __FILE__, __LINE__);
			return name;
		}
	}
	
	/*!
	 Returns the string "NULL" if failed, else user defined name of bond.
	 \param [in] index Internal index to locate and return the name associated.
	 */
	string System::bond_name (const int index) {
		string name = "NULL";
		char err_msg[MYERR_FLAG_SIZE];
		if (index >= bond_type_.size()) {
			sprintf(err_msg, "Index %d out of range", index);
			flag_error (err_msg, __FILE__, __LINE__);
			return name;
		} else {
			typedef std::map<std::string, int>::iterator it_type;
			for	(it_type iterator = bond_type_.begin(); iterator != bond_type_.end(); iterator++) {
				if (iterator->second == index) {
					return iterator->first;
				}
			}
			sprintf(err_msg, "Could not locate index %d in indexed bond types", index);
			flag_error (err_msg, __FILE__, __LINE__);
			return name;
		}
	}
	
	/*! Build interaction matrix for system
	  Right now, if two atoms are not bonded, then their interaction defaults to standard Lennard-Jones interaction
	*/
}		      

