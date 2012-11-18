/**
 MD Atom Information
 \author Nathan A. Mahynski
**/

#ifndef ATOM_H_
#define ATOM_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

//! Namespace containing pertinent Atom information
namespace atom {
	const int MAX_BONDS = 2;		//!< Maximum number of bonds any atom may have
	
	//! Atom class is defined so as to be easy to pass with MPI 
	class Atom {
		double pos[3];				//!< Cartesian coordinates
		double vel[3];				//!< Cartesian velocities
		double force[3];			//!< Cartesian force, (fx, fy, fz)
		double mass;				//!< Atomic mass (in reduced units)
		double diam;				//!< Atomic diameter (in reduced units)
		int type;					//!< Internally indexed type of this atom
		int	sys_index;				//!< Global atom index, i.e. unique in the system
		int nbonds;					//!< Number of bonds the atom has
		int bonds[MAX_BONDS];		//!< Global indices of bonds, fixed size for MPI passing
	};
	
	// Can define addition "atom-relevant" functions, classes, etc. below...
	
}

#endif