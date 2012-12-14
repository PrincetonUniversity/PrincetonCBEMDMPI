#include "CBEMD.h"
#include "mpi.h"
#include <limits>

//! Driver for mpi version
/*!
 \param \*argv[] ./cbemd nsteps dt xml_file bond_file
 */
int main (int argc, char *argv[]) {
	int rank, nprocs, rc;
	long int nsteps;
	double dt;
	System mysys;
	
	if (argc != 5) {
		fprintf(stderr, "syntax: ./cbemd nsteps dt xml_file energy_file\n");
		return ILLEGAL_VALUE;
	}
	
	// read nsteps and dt
	char* str_ptr;
	//	nsteps = atoi(argv[1]);
	nsteps = strtol(argv[1], &str_ptr, 10);
	if (str_ptr == argv[1]) {
	    fprintf(stderr, "No value was found for the number of steps. \n");
	    return ILLEGAL_VALUE;
	}
	if (*str_ptr != '\0') {
	    fprintf(stdout, "Number of steps taken as %ld ,", nsteps);
	    fprintf(stdout, " Ignored characters : %s\n", str_ptr);
	}

	dt = atof(argv[2]);
	if (nsteps < 0) {
		fprintf (stderr, "Error nsteps < 0\n");
		return ILLEGAL_VALUE;
	}
	if (dt < 0.0) {
		fprintf (stderr, "Error dt < 0.0\n");
		return ILLEGAL_VALUE;
	}
	
	// in the future "initialize_from_files" will take care of this line too
	Integrator *myint = new Verlet (dt);

	// set up MPI
	rc = MPI_Init(&argc, &argv);
	if (rc != MPI_SUCCESS) {
		printf ("Error starting MPI. Terminating.\n");
		MPI_Abort(MPI_COMM_WORLD, rc);
		return MPI_FAIL;
	}
	
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	// create and bind derived data type for atom
	create_MPI_ATOM();
	
	// initialize coordinates
	initialize_from_files (argv[3], argv[4], nprocs, rank, &mysys);
						 
	// initialize bonds/ppot
	//read_bond_file (argv[4], nprocs, rank, &mysys, myint);

	// set time step
	myint->set_dt(dt);
	
	// run
	run (&mysys, myint, nsteps);
	
	// print
	print_xml("outfile.xml", nprocs, rank, &mysys);
			  
	// finish
	delete_MPI_atom();
	MPI_Finalize();
	return SAFE_EXIT;
}
