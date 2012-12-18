/**
 Global variables
 \author Nathan A. Mahynski
 **/

#ifndef GLOBAL_H_
#define GLOBAL_H_

enum NAME_SIZES {ATOM_NAME_LENGTH = 100, BOND_NAME_LENGTH = 200, MYERR_FLAG_SIZE = 1000};
enum RETURN_FLAGS {SAFE_EXIT = 0, BAD_MEM = 1, ILLEGAL_VALUE = 2, MPI_FAIL = 3, INTEGRATE_FAIL = 4, FILE_ERROR = 5, BAD_EXIT = 6};				//!< Values that are returned if a failure condition is met
const double WCA_CUTOFF = pow(2.0, 1/6.0);
extern MPI_Datatype MPI_ATOM;

#endif
