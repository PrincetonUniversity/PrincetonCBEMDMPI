/**
 Miscellaneous Routines
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

//! Namespace for miscellaneous functions and tools
namespace misc {
	
	// George: define this at the top since this is just a header file 12/1/2012
	//const int ERR_FLAG_SIZE = 1000; //!< The maximum size allowed in an error or warning flag 
	
	//! Report an error message
	/*!
	 Error messages are piped to stderr not stdout.
	 \param [in] \*msg Character string to print out.
	 \param [in] \*file __FILE__ this function is called from.
	 \param [in] line __LINE__ this function is called from.
	 */
	void flag_error (const char *msg, const char *file, const int line);
    
	
	//! Report a notification
	/*!
	 Notification messages are piped to stderr not stdout.
	 \param [in] \*msg Character string to print out.
	 \param [in] \*file __FILE__ this function is called from.
	 \param [in] line __LINE__ this function is called from.
	 */
	void flag_notify (const char *msg, const char *file, const int line);
	
	//! Safely open a file
	/*!
	 Tries to open a file, if fails it returns a NULL pointer and alerts the user with an error message.  If it suceeds it returns the file pointer.
	 \param [in] \*filename Character name of file to open.
	 \param [in] \*opt File option ("r","w","rw+",etc.).
	 */
	FILE *mfopen(const char *filename, const char *opt);

	
	//! Returns the equivalent cartesian coordinates back in the simulation box assuming periodic boundaries.
	/*!
	 Recall that the simulation box is defined such that, regardless of the input file, it is normalized to 
	 have a corner at (0,0,0).  If this routine fails, it returns an empty vector (size = 0).
	 \param [in] coords Vector of cartesian coordinates.
	 \param [in] box Vector of box dimensions (L_x, L_y, L_z).
	 */
	vector <double> pbc (const vector <double> coords, const vector <double> box);
	
	
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


