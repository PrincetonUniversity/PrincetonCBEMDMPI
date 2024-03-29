/**
 \file initialize.cpp
 \brief Initialization routines for System object, as well as finalization of MPI.
 \author Nathan A. Mahynski
 **/

#include "initialize.h"

/*!
 Parse an XML and energy file to obtain atom and interaction information. Returns 0 if successful, non-zero if failure. Operates in
 a "cascade" between ranks so that each processor (if MPI is used) opens and initializes from the
 input file in order. 
 \param [in] xml_filename Name of coordinate file to open and read.
 \param [in] energy_filename Name of file containing bonds, pair potential parameters, etc.
 \param [in,out] \*sys Pointer to System object to store this information in.
 */
int initialize_from_files (const string xml_filename, const string energy_filename, System *sys) {
	int isignal, check, check2, check_sum, nprocs, rank;
	MPI_Status Stat;

	MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	MPI_Barrier (MPI_COMM_WORLD);
	
	// Cascade of reading statements
	if (rank == 0) {
		// read first, then signal next
		check = read_xml (xml_filename, sys);
		check2 = read_interactions (energy_filename, sys);
		check_sum = check+check2;
		
		if (nprocs > 1) {
			MPI_Send (&check_sum, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
		}
	} else {
		// wait mode for signal from rank-1
		MPI_Recv (&isignal, 1, MPI_INT, rank-1, rank-1, MPI_COMM_WORLD, &Stat);
		if (isignal == 0) {
			// recieved signal, read
			check = read_xml (xml_filename, sys);
			check2 = read_interactions (energy_filename, sys);
			check_sum = check+check2;
			
			// signal finished to next processor in line
			if (rank+1 < nprocs) {
				MPI_Send (&check_sum, 1, MPI_INT, rank+1, rank, MPI_COMM_WORLD);
			}
		} else {
			check_sum = isignal;
			if (rank+1 < nprocs) {
				MPI_Send (&isignal, 1, MPI_INT, rank+1, rank, MPI_COMM_WORLD);
			}
		}
	}
	
	MPI_Barrier (MPI_COMM_WORLD);
	return check_sum;
}

/*!
 Call MPI_Init and start the MPI ensuring it began successfully.  Returns SAFE_EXIT if successful.
 \param [in] argc Number of arguments in \*argv[].
 \param [in] \*argv[] Array of character arguments.
 */
int start_mpi (int argc, char *argv[]) {
	// set up MPI
	int rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		fprintf (stderr, "Error starting MPI. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		return MPI_FAIL;
	}
	
	// create and bind derived data type for atom
	create_MPI_ATOM();
	return SAFE_EXIT;
}

/*!
 Finalize MPI and clean up derived data types.
 */
int end_mpi () {
	delete_MPI_atom();
	MPI_Finalize();
	return SAFE_EXIT;
}

/*!
 Abort the MPI as a result of a failure, freeing MPI_ATOM as necessary.
 */
int abort_mpi () {
    delete_MPI_atom();
	MPI_Abort(MPI_COMM_WORLD, MPI_FAIL);
	return BAD_EXIT;
}
