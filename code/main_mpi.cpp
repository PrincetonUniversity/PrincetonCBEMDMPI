#include "CBEMD.h"
#include "mpi.h"

//! Driver for mpi version
/*!
 \param \*argv[] ./cbemd nsteps dt xml_file bond_file
 */
int main (int argc, char *argv[]) {
	int rank, nprocs, nsteps, rc;
	double dt;
	System mysys;
	Integrator *myint;
	
	if (argc != 5) {
		fprintf(stderr, "syntax: ./cbemd nsteps dt xml_file bond_file");
		return ILLEGAL_VALUE;
	}
	
	// read nsteps and dt
	nsteps = atoi(argv[1]);
	dt = atof(argv[2]);
	if (nsteps < 0) {
		fprintf (stderr, "Error nsteps < 0\n");
		return ILLEGAL_VALUE;
	}
	if (dt < 0.0) {
		fprintf (stderr, "Error dt < 0.0\n");
		return ILLEGAL_VALUE;
	}
	
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
	//initialize_from_xml (argv[3], nprocs, rank, &mysys);
						 
	// initialize bonds/ppot
	//read_bond_file (argv[4], nprocs, rank, &mysys, myint);

	// set time step
	myint->set_dt(dt);
	
	// run
	run (&mysys, myint, nsteps);
	
	// print?
	//print_xml("outfile.xml", nprocs, rank, &mysys);
			  
	// finish
	delete_MPI_atom();
	MPI_Finalize();
	return SAFE_EXIT;
}