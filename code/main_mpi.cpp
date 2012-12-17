#include "CBEMD.h"

//! Driver for mpi version
/*!
 \param \*argv[] ./cbemd nsteps dt xml_file bond_file
 */
int main (int argc, char *argv[]) {
	int rank, nprocs, rc, check;
	long int nsteps;
	double dt;
	System mysys;
	
	if (argc != 5) {
		fprintf(stderr, "syntax: ./cbemd nsteps dt xml_file energy_file\n");
		return ILLEGAL_VALUE;
	}
	
	// read nsteps and dt
	char* str_ptr;
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
	
	// set integrator
	Integrator *myint = new Verlet (dt);
	myint->set_dt(dt);
	
	check = start_mpi (argc, argv);
	if (check != 0) {
		goto finalize;
	}
	
	// initialize coordinates
	check = initialize_from_files (argv[3], argv[4], &mysys);
	if (check != 0) {
		goto finalize;
	}
	
	/*
	
	// run
	run (&mysys, myint, nsteps);
	
	// print
	print_xml("outfile.xml", nprocs, rank, &mysys);*/
			  
	// finish
finalize:
	end_mpi();
	return SAFE_EXIT;
}
