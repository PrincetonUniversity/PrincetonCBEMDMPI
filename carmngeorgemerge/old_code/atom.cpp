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
	const int num_fields = 9;

	MPI_Datatype type[num_fields] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_UB};
	int blocklen[num_fields] = {NDIM, NDIM, NDIM, NDIM, 1, 1, 1, 1, 1};
	Atom atom[2];
	MPI_Aint disp[num_fields];
	MPI_Aint start_address, address;
 
	MPI_Get_address(&(atom[0]), &start_address);
	
	MPI_Get_address(&(atom[0].pos[0]), &address);
	disp[0] = address - start_address;
	
	MPI_Get_address(&(atom[0].prev_pos[0]), &address);
	disp[1] = address - start_address;
	
	MPI_Get_address(&(atom[0].vel[0]), &address);
	disp[2] = address - start_address;
	
	MPI_Get_address(&(atom[0].force[0]), &address);
	disp[3] = address - start_address;

	MPI_Get_address(&(atom[0].mass), &address);
	disp[4] = address - start_address;
	
	MPI_Get_address(&(atom[0].diam), &address);
	disp[5] = address - start_address;
	
	MPI_Get_address(&(atom[0].type), &address);
	disp[6] = address - start_address;
	
	MPI_Get_address(&(atom[0].sys_index), &address);
	disp[7] = address - start_address;
	
	MPI_Get_address(&(atom[1]), &address);
	disp[8] = address - start_address;

	MPI_Type_create_struct(num_fields, blocklen, disp, type, &MPI_ATOM);
	MPI_Type_commit(&MPI_ATOM); 

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
