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
#include "interaction.h"
#include <list>
#include <algorithm>
#include <unordered_map>

using namespace std;
using namespace atom;
using namespace bond;
using namespace misc;

//! Namespace containing all system variable and quantities
namespace sim_system {
  	
  class System {
  public:
	  System();
	  ~System();
	  int set_box (const vector<double> new_box) {box_ = new_box;}	//!< Set the global system box size
	  vector <double> box() const {return box_;}				//!< Report system size
	  
	  void set_T (const double T) {Temp_ = T;}					//!< Set the system temperature
	  void set_P (const double P) {Press_ = P;}					//!< Set the system pressure
	  double T() const {return Temp_;}							//!< Report the temperature of the system
	  double P() const {return Press_;}							//!< Report the pressure of the system
	  
	  int natoms () const {return atoms_.size();}				//!< Return the number of atoms currently in this system (processor)
	  int add_atom_type (const string atom_name);				//!< Index an atom name
	  inline int atom_type (const string atom_name);			//!< Return the internal index associated with an atom name
	  inline string atom_name (const int index);				//!< Return the name associated with an index for atom type
    
	  int add_bond_type (const string bond_name);				//!< Index a bond name
	  inline int bond_type (const string bond_name);			//!< Return the internal index associated with a bond name
	  inline string bond_name (const int index);				//!< Return the name associated with an index for bond type
	  const pair <int, int> get_bond (const int nbond) {return bonded_[nbond];}	//!< Return a specific bonded pair indices
	  const int get_bond_type (const int nbond) {return bonded_type_[nbond];}	//!< Return the internal index of a bond
	  const int nbonds () {return bonded_.size();}				//!< Return the number of bonds in the system
	  inline void add_bond (const int atom1, const int atom2, const int type);
	
	  int* add_atoms (const int natoms, Atom *new_atoms);		//!< Add atom(s) to the system 
	  int delete_atoms (vector <int> indices);					//!< Pop atoms with local indices from local storage
	  Atom *get_atom (int index) {return &atoms_[index];}		//!< Get pointer to atom by local index
	

  private:
	  double Temp_;										//!< System temperature in reduced units (kT)
	  double Press_;									//!< System pressure in reduced units
	  double KE_, U_;									//!< Total internal kinetic energy and potential energy
	  vector <double> box_;								//!< Global system cartesian dimensions
	  vector <Atom> atoms_;								//!< Vector of Atoms in the system
	  vector <pair <int, int> > bonded_;				//!< Vector of global atom indices that are bonded (one record per bond)
	  vector <int> bonded_type_;						//!< Vector of bond type of each bond in bonded_ (in order)
	  unordered_map <string, int> atom_type_;			//!< Maps user specified name of atom type to internal index
	  unordered_map <string, int> bond_type_;			//!< Maps user specified name of bond type to internal index
	  unordered_map <string, int> ppot_type_;			//!< Maps user specified name of pair potential type to an internal index
	  unordered_map <int, int> glob_to_loc_id_;			//!< Maps global sys_index to the local index of atoms_ an atom is stored at on each processor; the opposite conversion can be done with lookup of Atom::sys_index
	  vector <int> on_proc_;							//!< Stores processor each atom (by global sys_index) "lives" on
	  int num_neighbor_procs;							//!< Number of neighboring processors this system must communicate with
	  int *neighbor_procs_;								//!< Array of ranks of neighboring processors this system must communicate with
	  vector <vector <Interaction> > interact_;			//!< Interaction matrix between atoms indexed by global id's (symetric)
  };
	
	
}

#endif

