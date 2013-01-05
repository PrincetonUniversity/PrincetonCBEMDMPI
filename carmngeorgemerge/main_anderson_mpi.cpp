/*!
 \brief Driver for MPI Version of CBEMD
 \file main_mpi.cpp
 \authors{Nathan A. Mahynski, Carmeline Dsilva, Arun L. Prabhu, George Khoury}
 **/

#include "CBEMD.h"

/*!
 \brief Main 
 \param \*argv[] ./andersen nsteps dt temperature xml_file bond_file output_file temperature
 */
int main (int argc, char *argv[]) {
	int check;
	long int nsteps;
	double dt;
	double temp;
	System mysys; //Declare system
	
	if (argc != 7) {
		fprintf(stderr, "syntax: ./andersen nsteps dt xml_file energy_file output_file temperature \n");
		return ILLEGAL_VALUE;
	}
	
	// read nsteps and dt
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

	// read temperature
	temp = atof(argv[6]);

	// setup integrator
	Integrator *myint = new Andersen (dt,temp);
	myint->set_dt(dt);  // Set timestep of integrator
	myint->set_temp(temp); // set temperature of integrator 

	check = start_mpi (argc, argv);
	if (check != 0) {
		goto finalize;
	}
	
	// initialize coordinates
	check = initialize_from_files (argv[3], argv[4], &mysys);
	if (check != 0) {
		goto finalize;
	}
	
	// run
	check = run (&mysys, myint, nsteps, argv[5]);
	if (check != 0) {
		goto finalize;
	}
	
	// print
	print_xml("outfile.xml", &mysys);
		  
	// finish
	end_mpi();
	return SAFE_EXIT;
	
finalize:
	abort_mpi();
	return BAD_EXIT;
}
