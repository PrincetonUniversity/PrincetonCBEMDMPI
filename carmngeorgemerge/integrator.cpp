/*!
   \brief MD Integrator(s) Information
   \authors{George Khoury, Carmeline Dsilva, Nathan A. Mahynski}
   \file integrator.cpp
**/

#include "integrator.h"

/*!
 \param [in] deltat Incremental timestep
 \param [in] temp Set point temperature
 */
Andersen::Andersen(double deltat, double temp, double nu){
	timestep_ = 0;
	dt_ = deltat;
	dt2_ = dt_*dt_;
	temp_ = temp;
	nu_ = nu;
}
	
/*!
 Step forward one timestep with the Andersen thermostat.  This reports the instantaneous temperature after each step.
 \param [in] \*sys Pointer to System to integrate.
 */
int Andersen::step (System *sys) {
	double prev_prev_pos;
	vector <double> box = sys->box();
	int check = 0, nprocs, rank;
	MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);

	double totalmass = 0;  
	for (int i = 0; i < sys->natoms(); ++i) {
		totalmass += sys->get_atom(i)->mass;
	}

	// Set the seed for the random number generator
	typedef std::tr1::linear_congruential<int, 16807, 0, (int)((1U << 31) -1 ) > Myceng;
	Myceng eng;

	// Step 1 is to use "velocity verlet" to integrate the positions
	for (int i = 0; i < sys->natoms(); ++i) {
		for (int j = 0; j < NDIM; ++j) {
			sys->get_atom(i)->prev_pos[j] = sys->get_atom(i)->pos[j];
			sys->get_atom(i)->pos[j] += sys->get_atom(i)->vel[j] * dt_ + 0.5 * sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
			sys->get_atom(i)->vel[j] += 0.5* dt_*sys->get_atom(i)->force[j] / sys->get_atom(i)->mass;
		}
	}

	// Step 2 is to find the forces after updating the positions and velocities
	check = force_calc(sys);
	double tempa = 0;
	for (int i = 0; i < sys->natoms(); ++i) {
		for (int j = 0; j < NDIM; ++j) {
			prev_prev_pos = sys->get_atom(i)->prev_pos[j];
			sys->get_atom(i)->prev_pos[j] = sys->get_atom(i)->pos[j];
			sys->get_atom(i)->vel[j] = sys->get_atom(i)->vel[j] + 0.5 * dt_ * sys->get_atom(i)->force[j] / sys->get_atom(i)->mass;
			tempa += sys->get_atom(i)->mass*sys->get_atom(i)->vel[j]*sys->get_atom(i)->vel[j];
		}
	}
			
	// We need the velocities from all the atoms from all procs together to get the instantaneous temperature
	double totaltempa;
	MPI_Allreduce (&tempa, &totaltempa,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	totaltempa = totaltempa / ( 3.0 * (sys->global_atom_types.size()));
	if (rank == 0) {
		cout << "Instantaneous Temp = "  << totaltempa << endl;
	}
			
	// Step 3 is to reset a certain number of velocities according the gaussian distribution
	double sig = sqrt(temp_);
	std::tr1::normal_distribution<double> distribution(0.0,sig);
	double rannum;
	for (int i = 0; i < sys->natoms(); ++i) {
		rannum = unifRand();
		if (rannum < nu_*dt_) {
			for (int j = 0; j < NDIM; ++j) {
				sys->get_atom(i)->vel[j] = distribution(eng);
			}
		}
	}
	timestep_++;
	return SAFE_EXIT;
}

/*!
 \param [in] deltat Incremental timestep
 */
Verlet::Verlet (double deltat) {
	timestep_ = 0;
	dt_ = deltat;
	dt2_ = dt_ * dt_;
	temp_ = -1.0;
}

/*! 
 \param [in,out] \*sys Pointer to System to make an integration step in.
*/
int Verlet::step (System *sys) {
	double prev_prev_pos;
	vector <double> box = sys->box();
	
	// On the first step, use euler-like step
	if (timestep_ == 0) {
		for (int i = 0; i < sys->natoms(); ++i) {
			for (int j = 0; j < NDIM; ++j) {
				sys->get_atom(i)->prev_pos[j] = sys->get_atom(i)->pos[j];
				sys->get_atom(i)->pos[j] += sys->get_atom(i)->vel[j] * dt_ + 0.5 * sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
				sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - sys->get_atom(i)->prev_pos[j]) / dt_;
			}
		}
	} else {
		for (int i = 0; i < sys->natoms(); ++i) {
			for (int j = 0; j < NDIM; ++j) {
				prev_prev_pos = sys->get_atom(i)->prev_pos[j];
				sys->get_atom(i)->prev_pos[j] = sys->get_atom(i)->pos[j];
				sys->get_atom(i)->pos[j] = 2.0 *  sys->get_atom(i)->prev_pos[j] - prev_prev_pos + sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
				sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - sys->get_atom(i)->prev_pos[j]) / dt_;
			}
		}
	}
	timestep_++;
	return SAFE_EXIT;
}
	
/*!
 Function returns SAFE_EXIT if successful, error flag otherwise.
 \param [in,out] \*sys Pointer to System to integrate
 \param [in] \*integrator Pointer to Integrator to use
 \param [in] timesteps Number of timesteps to integrate over
 \param [in] outfile Name of file to print animation to
*/
int run (System *sys, Integrator *integrator, const int timesteps, const string outfile) {
	const bool wrap_pos = true;
	int check = 0, nprocs, rank;
	char err_msg[MYERR_FLAG_SIZE];
  
	MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
		
	/* Before starting, we need to check that all requisite variables are set */
	// Check interactions have been properly set
	if (sys->interact.size() == 0) {
		sprintf(err_msg, "Interactions have not been set on rank %d", rank);
		flag_error (err_msg, __FILE__, __LINE__);
		return ILLEGAL_VALUE;
	} else {
		for (unsigned int i = 0; i < sys->interact.size(); ++i) {
			if (sys->interact[i].size() != sys->interact.size()) {
				sprintf(err_msg, "Interactions have not been fully set on rank %d", rank);
				flag_error (err_msg, __FILE__, __LINE__);
				return ILLEGAL_VALUE;
			}
			for (unsigned int j = 0; j < sys->interact.size(); ++j) {
				if (sys->interact[i][j].check_force_energy_function() == NULL && i != j) {
					sprintf(err_msg, "Interactions have not been fully set on rank %d", rank);
					flag_error (err_msg, __FILE__, __LINE__);
					return ILLEGAL_VALUE;
				}
			}
		}
	}
		
	// Check their are some atoms in the system
	if (sys->global_atom_types.size() < 1) {
		sprintf(err_msg, "No atoms found in system");
		flag_error (err_msg, __FILE__, __LINE__);
		return ILLEGAL_VALUE;
	}
		
	// Check integrator has been set
	if (integrator == NULL) {
		sprintf(err_msg, "No integrator set on rank %d", rank);
		flag_error (err_msg, __FILE__, __LINE__);
		return ILLEGAL_VALUE;
	}
		
	// Check timesteps, dt > 0
	if (integrator->dt() <= 0.0) {
		sprintf(err_msg, "dt < 0 on rank %d", rank);
		flag_error (err_msg, __FILE__, __LINE__);
		return ILLEGAL_VALUE;
	}
	if (timesteps < 0) {
		sprintf(err_msg, "Timsteps < 0 on rank %d", rank);
		flag_error (err_msg, __FILE__, __LINE__);
		return ILLEGAL_VALUE;
	}
	
	// Estimate how often to report the progress of the simulation
	int print_step;
	if (timesteps < 100) {
		print_step = 1;
	} else {
		print_step = (int) floor(timesteps/100.0);
	}

	// Execute loop
	for (int i = 0; i < timesteps; ++i) {
		// All steps must happen at the same time
		MPI_Barrier(MPI_COMM_WORLD);
		
		// Move atoms that have left the simulation box
		if (nprocs > 1) {
			check = send_atoms(sys);
			if (check != 0) {
				sprintf(err_msg, "Error encountered during sending atoms after step %d", i+1);
				flag_error (err_msg, __FILE__, __LINE__);
				return check;
			}
		}

		// Clear out the forces on each atom which is NECESSARY before each new step
		for (int j = 0; j < sys->natoms(); ++j) {
			for (int k = 0; k < NDIM; ++k) {
				sys->get_atom(j)->force[k] = 0.0;
			}
		}

		// Calculate forces
		if (i == 0) {
			// on first step, all integrators need the initial forces
			check = force_calc(sys);
		} else if (integrator->getTemp() < 0) {
			// If NVE is being used get force, otherwise this is done in integrator.step()
			check = force_calc(sys);
		} 
		
		if (check != 0) {
			sprintf(err_msg, "Error encountered during force calc after step %d", i+1);
			flag_error (err_msg, __FILE__, __LINE__);
			return check;
		}

		// Create animation
		write_xyz (outfile, sys, i, wrap_pos);

		// Step forward
		check = integrator->step(sys);
		if (check != 0) {
			sprintf(err_msg, "Error encountered during integration after step %d", i+1);
			flag_error (err_msg, __FILE__, __LINE__);
			return check;
		}

		// Report progress
		if (rank == 0) {
			if (i%print_step == 0) {
				sprintf(err_msg, "Finished %d of %d total steps", i, timesteps);
				flag_notify (err_msg, __FILE__, __LINE__);
			}
		}
	}
	
	// Report the final positions
	write_xyz (outfile, sys, timesteps, wrap_pos);
		
	return SAFE_EXIT;
}

