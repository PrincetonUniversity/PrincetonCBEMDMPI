/*!
 \file misc.h
 \brief Header for Miscellaneous Routines
 \author Nathan A. Mahynski
 **/

#ifndef MISC_H_
#define MISC_H_

#include "common.h"
#include <map>
#include <assert.h>
#include "atom.h"
#include "global.h"

using namespace std;
using namespace atom;

/* \namespace misc
  Namespace for miscellaneous functions and tools
*/
namespace misc {
		
	//! Report an error message
	void flag_error (const char *msg, const char *file, const int line);
    
	
	//! Report a notification
	void flag_notify (const char *msg, const char *file, const int line);
	
	//! Safely open a file
	FILE *mfopen(const char *filename, const char *opt);

	
	//! Returns the equivalent cartesian coordinates back in the simulation box assuming periodic boundaries.
	vector <double> pbc (const vector <double> coords, const vector <double> box);
	
	//! Returns the equivalent cartesian coordinates back in the simulation box assuming periodic boundaries.
	vector <double> pbc (const double *coords, const vector <double> box);
	
	//! Return the square of the minimum image distance between 2 coordinate vectors
	/*!
	 \param [in] coords1 Vector of cartesian coordinates of one atom
	 \param [in] coords2 Vector of cartesian coordinates of the other atom
	 \param [in] box Vector of cartesian coordinates of the box
	 */
	double min_image_dist2 (const vector <double> coords1, const vector <double> coords2, const vector <double> box); 
	
	//! Returns the square of the minimum image distance between 2 atoms, also returns the minimum image distance vector xyz that points from atom1 to atom2
	/*!
	 \param [in] \*a1 Pointer to one atom
	 \param [in] \*a2 Pointer to the other atom
	 \param [in] box Vector of cartesian coordinates of the box
	 \param [in,out] \*xyz Array of xyz displacements to be returned to the user (length 3)
	 */
	double min_image_dist2 (const Atom *a1, const Atom *a2, const vector <double> *box, double *xyz);
	
}


#endif


