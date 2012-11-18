/**
 Read HOOMD XML file format into System object
 \author Nathan A. Mahynski
 **/

#ifndef READ_XML_H_
#define READ_XML_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <map>
#include "system.h"
#include "misc.h"
#include "atom.h"

using namespace std;
using namespace system;
using namespace misc;
using namespace atom;

/*!
 Parse an XML file to obtain atom information. Returns 0 if successful, -1 if failure.
 \param [in] filename Name of file to open and read.
 */
int read_xml (const string filename, const int nprocs) {
	FILE *input = mfopen(filename, "r");
	
	
}

/*!
 \param [in] \*filename Name of file to initialize from
 \param [in[ rank Rank of processor this is.  Rank 0 reads, other wait to recieve information
 \param [in,out] \*sys System object to store information this processor is responsible for
 */
int initialize (const char *filename, int *rank, int *nprocs, System *sys) {
	int argc, rc;
	char *argv[];
	
	// set up MPI
	rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		return MPI_FAIL;
	}
	
	MPI_Comm_size(MPI_COMM_WORLD,nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,rank);
	
	// Create MPI_ATOM datatype
	create_MPI_ATOM();
	
	// check that interactions are such that ONLY neighboring procs need to interact
	
	// handle domain decomp --> look at MPI_graph?
	
	if (rank == 0) {
		int check = read_xml(	...	);
	} else {
		int check = worker_recv_sys(	...	);
	}
	
	return 0;
}

int finalize () {
	// free atom type after we are done running
	delete_MPI_atom();
	MPI_Finalize();
	return flag;
}
#endif