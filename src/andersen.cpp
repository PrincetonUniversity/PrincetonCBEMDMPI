/*!
 \brief Driver for MPI Version of CBEMD with Andersen Thermostat
 \file andersen.cpp
 \authors{Nathan A. Mahynski, Carmeline Dsilva, Arun L. Prabhu, George Khoury, Frank Ricci, Jun Park}
 **/

#include "CBEMD.h"

/*!
 \param \*argv[] ./andersen nsteps dt xml_file energy_file bond_file animation_file temperature nu
 */
int main (int argc, char *argv[]) {
	int check;
	long int nsteps;
	double dt;
	System mysys; //Declare system
	
	if (argc != 8) {
		fprintf(stderr, "syntax: ./andersen nsteps dt xml_file energy_file animation_file temperature nu \n");
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

	// Read temperature
	double temp = atof(argv[6]);
	
	// Read nu
	double nu = atof(argv[7]);

	// Setup integrator
	Integrator *myint = new Andersen (dt,temp,nu);
	myint->set_dt(dt);  
	myint->set_temp(temp); 

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
	
	// Print
	print_xml("outfile.xml", &mysys);
		  
	// Finish if successful
	end_mpi();
	return SAFE_EXIT;
	
	// Finish if error was reported
finalize:
	abort_mpi();
	return BAD_EXIT;
}
