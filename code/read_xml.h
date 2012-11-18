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

using namespace std;
using namespace system;
using namespace misc;

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
int initialize (const char *filename, const int rank, System *sys) {
	// check that interactions are such that ONLY neighboring procs need to interact
	
	// handle domain decomp --> look at MPI_graph?
	
	if (rank == 0) {
		int check = read_xml(	...	);
	} else {
		int check = worker_recv_sys(	...	);
	}
}

#endif