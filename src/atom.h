/*!
 \file atom.h
 \brief MD Atom Information
 \author Nathan A. Mahynski
**/

#ifndef ATOM_H_
#define ATOM_H_

#include "mpi.h"
#include "global.h"

//! Atom class is defined as a struct so as to be easy to pass with MPI
typedef struct {
	double pos[NDIM];			//!< Cartesian coordinates
	double prev_pos[NDIM];      //!< Cartesian coordinates for the previous position of the atom (needed for integrator)
	double vel[NDIM];			//!< Cartesian velocities (vx, vy, vz)
	double force[NDIM];			//!< Cartesian force, (fx, fy, fz)
	double mass;				//!< Atomic mass (in reduced units)
	double diam;				//!< Atomic diameter (in reduced units)
	int type;					//!< Internally indexed type of this atom
	int	sys_index;				//!< Global atom index, i.e. unique in the system
} Atom;

//! Creates the MPI_Atom class so it can be passed with MPI
void create_MPI_ATOM ();
	
//! Free the MPI type at the end of the program
void delete_MPI_atom();

#endif
