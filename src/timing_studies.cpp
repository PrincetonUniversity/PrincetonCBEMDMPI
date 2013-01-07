/*!
 \brief Driver for MPI Version of CBEMD
 \file main_mpi.cpp
 \authors{Nathan A. Mahynski, Carmeline Dsilva, Arun L. Prabhu, George Khoury}
 **/

#include "CBEMD.h"

int main (int argc, char *argv[]) {
	int check;
	long int nsteps = atoi(argv[1]);
	double dt = 0.003;
	System mysys; //Declare system
	
	string xml_file = argv[2];
	string energy_file = argv[3];
	string output_file = argv[4];

	// setup integrator
	Integrator *myint = new Verlet (dt);
	myint->set_dt(dt);  // Set timestep of integrator

	check = start_mpi (argc, argv);
	if (check != 0) {
		goto finalize;
	}
	
	// initialize coordinates
	check = initialize_from_files (xml_file, energy_file, &mysys);
	if (check != 0) {
		goto finalize;
	}
	
	// run
	check = run (&mysys, myint, nsteps, output_file);
	if (check != 0) {
		goto finalize;
	}
	
	// finish
	end_mpi();
	return SAFE_EXIT;
	
finalize:
	abort_mpi();
	return BAD_EXIT;
}