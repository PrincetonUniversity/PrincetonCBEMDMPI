/**
 \file global.h
 \brief Global variables
 \authors{Nathan A. Mahynski, George Khoury}
 **/

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <exception>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>

//! Set maximum size for user-defined names of things
enum NAME_SIZES {ATOM_NAME_LENGTH = 100, BOND_NAME_LENGTH = 200, MYERR_FLAG_SIZE = 1000};

//! Set error flags returned if a failure condition is met
enum RETURN_FLAGS {SAFE_EXIT = 0, BAD_MEM = 1, ILLEGAL_VALUE = 2, MPI_FAIL = 3, INTEGRATE_FAIL = 4, FILE_ERROR = 5, BAD_EXIT = 6};				

//! The Weeks-Chandler-Andersen cutoff length is a fixed number that is used in numerous places commonly in MD so it is precomputed here
const double WCA_CUTOFF = pow(2.0, 1/6.0);

//! MPI_ATOM is made visible to all routines
extern MPI_Datatype MPI_ATOM;

//! Number of dimensions in the system (3D), but future work may include 2D as well.  This leaves the door open for this.
const int NDIM = 3;  

//! Nearest neighbors for 3D Domain decomposition
const int NNEIGHBORS = 26;

//! dimension along which to do the domain decomposition
const int PARALLELDIM = 0;

#endif
