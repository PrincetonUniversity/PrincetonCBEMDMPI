/*! 
 \file misc.cpp
 \brief Source code for miscellaneous routines
 \author Nathan A. Mahynski
**/

#include "misc.h"
        
/*!
 Error messages are sent to stderr.
 \param [in] \*msg Character string to print out.
 \param [in] \*file __FILE__ this function is called from.
 \param [in] line __LINE__ this function is called from.
 */
void flag_error (const char *msg, const char *file, const int line) {			
	fprintf(stderr, "*** Error :: %s :: [FILE: %s, LINE %d] ***\n", msg, file, line);
}

/*!
 Notification messages are sent to stderr.
 \param [in] \*msg Character string to print out.
 \param [in] \*file __FILE__ this function is called from.
 \param [in] line __LINE__ this function is called from.
 */
void flag_notify (const char *msg, const char *file, const int line) {			
	fprintf(stderr, "--- Note :: %s :: [FILE: %s, LINE %d] ---\n", msg, file, line);
}

/*!
 Tries to open a file, if it fails it returns a NULL pointer and alerts the user with an error message.  If it succeeds it returns the file pointer.
 \param [in] \*filename Character name of file to open.
 \param [in] \*opt File option ("r","w","rw+",etc.).
 */
FILE *mfopen(const char *filename, const char *opt) {
	FILE *fp1;
	char err_msg[MYERR_FLAG_SIZE]; 
	if (!(fp1 = fopen(filename, opt))) {
		sprintf(err_msg, "Could not open %s with %s permissions", filename, opt);
		flag_error(err_msg, __FILE__, __LINE__);
		return NULL;
	}
	return fp1;
}

/*!
 If this routine fails, it returns an empty vector (size = 0).
 \param [in] coords Vector of cartesian coordinates.
 \param [in] box Vector of box dimensions (L_x, L_y, L_z).
 */
vector <double> pbc (const vector <double> coords, const vector <double> box) {
	vector <double> in_box(NDIM, 0.0); 
	for (int i = 0; i < NDIM; ++i) {
	    in_box[i] = coords[i];
		// because box is defined with corner at (0,0,0) ceil/floor when used appropriately puts back in box with one operation
		in_box[i] = coords[i];
		if (coords[i] < 0.0) {
			in_box[i] = coords[i] + ceil(-coords[i]/box[i])*box[i];
		}
		if (coords[i] >= box[i]) {
			in_box[i] = coords[i] - floor(coords[i]/box[i])*box[i];
		}
	}
	
	return in_box;
}
	
/*!
 If this routine fails, it returns an empty vector (size = 0).
 \param [in] \*coords Array of cartesian coordinates.
 \param [in] box Vector of box dimensions (L_x, L_y, L_z).
*/
vector <double> pbc (const double *coords, const vector <double> box) {
	vector <double> in_box(NDIM, 0.0);
	for (int i = 0; i < NDIM; ++i) {
		in_box[i] = coords[i];
		if (coords[i] < 0.0) {
			in_box[i] = coords[i] + ceil(-coords[i]/box[i])*box[i];
		}
		if (coords[i] >= box[i]) {
			in_box[i] = coords[i] - floor(coords[i]/box[i])*box[i];
		}
	}
	return in_box;
}

/*!
 \param [in] coords1 Vector of cartesian coordinates of one atom
 \param [in] coords2 Vector of cartesian coordinates of the other atom
 \param [in] box Vector of cartesian coordinates of the box
 */
double min_image_dist2 (const vector <double> coords1, const vector <double> coords2, const vector <double> box) {
	double ans = 0.0, dist;
	for (int i = 0; i < NDIM; ++i) {
		dist = coords2[i] - coords1[i];
		dist -= round(dist/box[i])*box[i];
		ans += dist*dist;
	}
	return ans;
}

/*!
 \param [in] \*a1 Pointer to one atom
 \param [in] \*a2 Pointer to the other atom
 \param [in] box Vector of cartesian coordinates of the box
 \param [in,out] \*xyz Array of xyz displacements to be returned to the user (length 3)
 */
double min_image_dist2 (const Atom *a1, const Atom *a2, const vector <double> *box, double *xyz) {
	double ans = 0.0;
	for (int i = 0; i < NDIM; ++i) {
		xyz[i] = a2->pos[i] - a1->pos[i];
		xyz[i] -= round(xyz[i]/box->at(i))*box->at(i);
		ans += xyz[i]*xyz[i];
	}
	return ans;
}

double unifRand() {
	return rand() / double(RAND_MAX);
}

/*! 
 \param a Lower end point of the interval
 \param b Upper end of the interval
 */
double unifRand(double a, double b) {
	return (b-a)*unifRand() + a;
}
