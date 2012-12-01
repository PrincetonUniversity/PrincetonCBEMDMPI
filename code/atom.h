/**
 MD Atom Information
 \author Nathan A. Mahynski
**/

#ifndef ATOM_H_
#define ATOM_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"

//! Namespace containing pertinent Atom information
namespace atom {
	//! Atom class is defined so as to be easy to pass with MPI 
	struct Atom {
		double pos[3];				//!< Cartesian coordinates
		double vel[3];				//!< Cartesian velocities
		double force[3];			//!< Cartesian force, (fx, fy, fz)
		double mass;				//!< Atomic mass (in reduced units)
		double diam;				//!< Atomic diameter (in reduced units)
		int type;					//!< Internally indexed type of this atom
		int	sys_index;				//!< Global atom index, i.e. unique in the system
	};
	
	//! MPI version of Atom
	MPI_Datatype MPI_ATOM;
	
	//! Creates the MPI_Atom class so it can be passed over MPI
	/*
	 \sa initialize
	 */
	// (see example of use at https://computing.llnl.gov/tutorials/mpi/#Derived_Data_Types)
	void create_MPI_ATOM () {
		int num_types = 2;
		MPI_Datatype oldtypes[num_types];
		MPI_Aint offsets[num_types], extent;
		int blockcounts[num_types];
		
		MPI_Type_extent (MPI_DOUBLE, &extent);
		
		offsets[0] = 0;
		blockcounts[0] = 11;
		oldtypes[0] = MPI_DOUBLE;
		
		offsets[1] = (11*extent);
		blockcounts[1] = 2;
		oldtypes[1] = MPI_INT;
		
		MPI_Type_struct (num_types, blockcounts, offsets, oldtypes, &MPI_ATOM);
		MPI_Type_commit (&MPI_ATOM);
	}
	
	//! Free the MPI type at the end of the program
	/*! 
	 \sa finalize
	 */
	void delete_MPI_atom() {
		MPI_Type_free (&MPI_ATOM);
	}
}

#endif