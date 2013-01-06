/*!
 \brief Driver for MPI Version of CBEMD verlet
 \file verlet.cpp
 \authors{Nathan A. Mahynski, Carmeline Dsilva, Arun L. Prabhu, George Khoury, Frank Ricci, Jun Park}
 **/

#include "CBEMD.h"

/*!
 \param \*argv[] ./verlet nsteps dt xml_file energy_file animation_file
 */
int main (int argc, char *argv[]) {
	int check;
	long int nsteps;
	double dt;
	System mysys; //Declare system
	
	if (argc != 6) {
		fprintf(stderr, "syntax: ./verlet nsteps dt xml_file energy_file animation_file\n");
		return ILLEGAL_VALUE;
	}
	
	// Read nsteps and dt
	char* str_ptr;
	nsteps = strtol(argv[1], &str_ptr, 10); // Store number of steps
	dt = atof(argv[2]);
	if (str_ptr == argv[1]) {
	    fprintf(stderr, "No value was found for the number of steps. \n");
	    return ILLEGAL_VALUE;
	}
	if (*str_ptr != '\0') {
	    fprintf(stdout, "Number of steps taken as %ld ,", nsteps);
	    fprintf(stdout, " Ignored characters : %s\n", str_ptr);
	}

	// Setup integrator
	Integrator *myint = new Verlet (dt);
	myint->set_dt(dt);  

	check = start_mpi (argc, argv);
	if (check != 0) {
		goto finalize;
	}
	
	// Initialize coordinates
	check = initialize_from_files (argv[3], argv[4], &mysys);
	if (check != 0) {
		goto finalize;
	}
	
	// Run
	check = run (&mysys, myint, nsteps, argv[5]);
	if (check != 0) {
		goto finalize;
	}
	
	// Print final coordinates
	print_xml("outfile.xml", &mysys);
		  
	// Finish safely
	end_mpi();
	return SAFE_EXIT;
	
	// Finish if error was reported
finalize:
	abort_mpi();
	return BAD_EXIT;
}
