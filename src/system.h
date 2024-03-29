/**
 \file system.h
 \brief Header for MD System Information
 \author Nathan A. Mahynski
 **/

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <map>
#include "atom.h"
#include "misc.h"
#include "interaction.h"
#include <list>
#include <algorithm>
#include "global.h"

using namespace std;

class System {
public:
	System();
	~System();
	void set_box (const vector<double> new_box) {box_ = new_box;}	//!< Set the global system box size
	vector <double> box() const {return box_;}				//!< Report system size
		
	void set_T (const double T) {Temp_ = T;}				//!< Set the system temperature
	void set_P (const double P) {Press_ = P;}				//!< Set the system pressure
	double T() const {return Temp_;}						//!< Report the temperature of the system
	double P() const {return Press_;}						//!< Report the pressure of the system
		
	int total_atoms () const {return atoms_.size();}		//!< Return the number of atoms currently in this system (processor) including current ghosts
	int natoms () const {return num_atoms_;}				//!< Return the number of atoms this system (processor) is responsible for
	int add_atom_type (const string atom_name);				//!< Index an atom name
	int atom_type (const string atom_name);					//!< Return the internal index associated with an atom name
	string atom_name (const unsigned int index);			//!< Return the name associated with an index for atom type

	int add_bond_type (const string bond_name);				//!< Index a bond name
	int bond_type (const string bond_name);					//!< Return the internal index associated with a bond name
	string bond_name (const unsigned int index);								//!< Return the name associated with an index for bond type
	const pair <int, int> get_bond (const int nbond) {return bonded_[nbond];}	//!< Return a specific bonded pair indices
	int get_bond_type (const int nbond) {return bonded_type_[nbond];}			//!< Return the internal index of a bond
	int nbonds () {return bonded_.size();}					//!< Return the number of bonds in the system
	void add_bond (const int atom1, const int atom2, const int type);
		
	vector <int> add_atoms (const int natoms, Atom *new_atoms);		//!< Add atom(s) to the system with an array of atoms
	void add_ghost_atoms (const int natoms, Atom *new_atoms);		//!< Add ghost atom(s) to the system (does not update the number of atoms the processor is responsible for
	vector <int> add_atoms (vector <Atom> *new_atoms);				//!< Add atom(s) to the system with an vector of atoms
	int delete_atoms (vector <int> indices);				//!< Pop atoms with local indices from local storage
	Atom *get_atom (int index) {return &atoms_[index];}		//!< Get pointer to atom by local index
	Atom copy_atom (int index) {return atoms_[index];}		//!< Report a copy of an atom
	void set_rank (int rank) {rank_ = rank;}				//!< Record the rank this system corresponds to
	int rank () {return rank_;}								//!< Return the rank of the system
	void set_num_atoms (int size) {num_atoms_ = size;}		//!< Manually set the number of atoms in the system
	void clear_ghost_atoms ();								//!< Clear ghost atoms from system
		
	/* These are associated with 3D Domain Decomp */
	int gen_domain_info ();
	double proc_widths[NDIM];								//!< Width for domain decomposition
	vector<int> final_proc_breakup;							//!< Final domain decomposition
	int xyz_id[NDIM];
	double xyz_limits[NDIM][2];
	int send_table [NNEIGHBORS];
	vector< vector<Atom> > send_lists;
	int send_list_size[NNEIGHBORS], get_list_size[NNEIGHBORS];
	vector< vector<Atom> > get_lists;
		
	vector <vector <Interaction> > interact;				//!< Interaction matrix between atoms indexed by global id's (symetric)
	vector <string> global_atom_types;						//!< Keeps a record of every atom's type

	void set_max_rcut (const double max_rcut) {max_rcut_ = max_rcut;}	//!< Set the maximum cutoff radius of all interactions in the system
	double max_rcut () const {return max_rcut_;}						//!< Return the max cutoff radius
		
	void set_total_KE (const double ke) {KE_ = ke;}			//!< Set the global kinetic energy record
	void set_total_PE (const double pe) {U_ = pe;}			//!< Set the global potential energy record
	double KE () const {return KE_;}
	double U () const {return U_;}
		
private:
	int rank_;										//!< Rank of the processor this domain is on
	double Temp_;									//!< System temperature in reduced units (kT)
	double Press_;									//!< System pressure in reduced units
	double KE_, U_;									//!< Total internal kinetic energy and potential energy of the global system
	vector <Atom> atoms_;							//!< Vector of Atoms in the system
	vector <double> box_;							//!< Global system cartesian dimensions
	vector <double> masses_;						//!< Vector containing masses of each type of Atom in the system
	map <string, unsigned int> atom_type_;			//!< Maps user specified name of atom type to internal index
	map <string, unsigned int> bond_type_;			//!< Maps user specified name of bond type to internal index
	map <string, unsigned int> ppot_type_;			//!< Maps user specified name of pair potential type to an internal index
	map <int, int> glob_to_loc_id_;					//!< Maps global sys_index to the local index of atoms_ an atom is stored at on each processor; the opposite conversion can be done with lookup of Atom::sys_index
	vector < pair <int, int> > bonded_;				//!< Vector of bonded pairs
	vector <int> bonded_type_;						//!< Vector of types associated with each bond
	int num_atoms_;									//!< The number of atoms the processor is responsible for
	double max_rcut_;								//!< Max cutoff radius for all interactions used in the system
};

#endif
