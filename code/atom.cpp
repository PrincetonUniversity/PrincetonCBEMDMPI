/*!
 \file atom.cpp
 \brief Functions for handling MD Atom Information
\author Nathan A. Mahynski
**/


#include "atom.h"

//! \namespace atom
namespace atom {

/*!
 * \sa initialize
 * This function creates and commits to the system the MPI_ATOM derived data type.
 * (see example of use at https://computing.llnl.gov/tutorials/mpi/#Derived_Data_Types)
 */
void create_MPI_ATOM () {
	const int num_types = 2;
	MPI_Datatype oldtypes[num_types];
	MPI_Aint offsets[num_types], extent;
	int blockcounts[num_types];
	//	extern MPI_Datatype MPI_ATOM;

	MPI_Type_extent (MPI_DOUBLE, &extent);
	
	offsets[0] = 0;
	blockcounts[0] = 14;
	oldtypes[0] = MPI_DOUBLE;
	
	offsets[1] = (14*extent);
	blockcounts[1] = 2;
	oldtypes[1] = MPI_INT;
	
	MPI_Type_struct (num_types, blockcounts, offsets, oldtypes, &MPI_ATOM);
	MPI_Type_commit (&MPI_ATOM);
}

/*!
 * This function utilizes the complementary MPI_Type_free routine
 * to mark the MPI_ATOM type for deallocation
 * \sa finalize
 */
void delete_MPI_atom() {
	MPI_Type_free (&MPI_ATOM);
}

}
