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
#include "global.h"

//! Namespace containing pertinent Atom information
namespace atom {
	//! Atom class is defined so as to be easy to pass with MPI 
	typedef struct {
		double pos[3];				//!< Cartesian coordinates
		double vel[3];				//!< Cartesian velocities
		double force[3];			//!< Cartesian force, (fx, fy, fz)
		double mass;				//!< Atomic mass (in reduced units)
		double diam;				//!< Atomic diameter (in reduced units)
		int type;					//!< Internally indexed type of this atom
		int	sys_index;				//!< Global atom index, i.e. unique in the system
	} Atom;
	
	//! MPI version of Atom
	//extern MPI_Datatype MPI_ATOM; // George : is this global? Maybe we should put in driver. 11/30 11:51PM
	
	//! Creates the MPI_Atom class so it can be passed over MPI
	/*
	 \sa initialize
	 */
	// (see example of use at https://computing.llnl.gov/tutorials/mpi/#Derived_Data_Types)
	void create_MPI_ATOM ();
	
	//! Free the MPI type at the end of the program
	/*! 
	 \sa finalize
	 */
	void delete_MPI_atom();
}

#endif
