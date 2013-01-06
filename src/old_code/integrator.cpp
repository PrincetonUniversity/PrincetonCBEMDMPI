/**
   \brief MD Integrator(s) Information
   \authors{George Khoury, Carmeline Dsilva, Nathan A. Mahynski}
   \file integrator.cpp
**/

#include "integrator.h"
#include "force_calc.h"

//! \namespace integrator
//! Namespace for the base class for integrators
namespace integrator {

	//! NVE, Verlet
	Verlet::Verlet (double deltat) {
		timestep_ = 0;
		dt_ = deltat;
		dt2_ = dt_ * dt_;
	}

	/*! 
	  This function integrates.
	  Call it to make a step.
	  Will update atom positions and velocities
	  Maintains atom positions in the box if atoms cross boundaries.
	  \param [in,out] sys Pointer to System to make an integration step in
	*/
	int Verlet::step (System *sys) {
		double prev_prev_pos;
		vector <double> box = sys->box();
		if (timestep_ == 0) {
			for (int i = 0; i < sys->natoms(); ++i) {
				for (int j = 0; j < NDIM; ++j) {
					sys->get_atom(i)->prev_pos[j] = sys->get_atom(i)->pos[j];
					sys->get_atom(i)->pos[j] += sys->get_atom(i)->vel[j] * dt_ + 0.5 * sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
					// maintain atom position in the box
					// use direct iteration NOT ceil/floor here because positions should not move excessively unless you have bigger problems
					/*while (sys->get_atom(i)->pos[j] < 0.0) {
						sys->get_atom(i)->pos[j] += box[j];
					}
					while (sys->get_atom(i)->pos[j] > box[j]) {
						sys->get_atom(i)->pos[j] -= box[j];
					}*/
					/*if (sys->get_atom(i)->pos[j] < 0.0) {
						sys->get_atom(i)->pos[j] += ceil(-sys->get_atom(i)->pos[j]/box[j])*box[j];
					}
					if (sys->get_atom(i)->pos[j] >= box[j]) {
						sys->get_atom(i)->pos[j] -= floor(sys->get_atom(i)->pos[j]/box[j])*box[j];
						}*/
					sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - sys->get_atom(i)->prev_pos[j]) / dt_;
				}
			}
		}
		else {
			for (int i = 0; i < sys->natoms(); ++i) {
				for (int j = 0; j < NDIM; ++j) {
					prev_prev_pos = sys->get_atom(i)->prev_pos[j];
					sys->get_atom(i)->prev_pos[j] = sys->get_atom(i)->pos[j];
					sys->get_atom(i)->pos[j] = 2.0 *  sys->get_atom(i)->prev_pos[j] - prev_prev_pos + sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
					// maintain atom position in the box
					// use direct iteration NOT ceil/floor here because positions should not move excessively unless you have bigger problems
					/*while (sys->get_atom(i)->pos[j] < 0.0) {
						sys->get_atom(i)->pos[j] += box[j];
					}
					while (sys->get_atom(i)->pos[j] > box[j]) {
						sys->get_atom(i)->pos[j] -= box[j];
					}*/
					/*if (sys->get_atom(i)->pos[j] < 0.0) {
						sys->get_atom(i)->pos[j] += ceil(-sys->get_atom(i)->pos[j]/box[j])*box[j];
					}
					if (sys->get_atom(i)->pos[j] >= box[j]) {
						sys->get_atom(i)->pos[j] -= floor(sys->get_atom(i)->pos[j]/box[j])*box[j];
						}*/
					//sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - prev_prev_pos) / (2.0 * dt_);
					sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - sys->get_atom(i)->prev_pos[j]) / dt_;
				}
			}
		}
		timestep_++;
		return 0;
	}
	
	
	//!< Run (i.e. integrate) a system forward in time for a specified number of timesteps
	/*!
	  Function returns 0 if successful, -1 if it encountered an error.
	  \param [in] sys Pointer to System to integrate
	  \param [in] integrator Pointer to Integrator to use
	  \param [in] timesteps Number of timesteps to integrate over
	*/
	int run (System *sys, Integrator *integrator, const int timesteps, const string outfile) {
		const bool wrap_pos = true;

		// The way this function is written it can be easily interpreted by SWIG with python!
		int check = 0, nprocs, rank;
		char err_msg[MYERR_FLAG_SIZE];
  
		MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
		MPI_Comm_rank (MPI_COMM_WORLD, &rank);
		
		/* before starting, need to check that all requisite variables are set */
		// check interactions have been properly set
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
		
		// check their are some atoms in the system
		if (sys->global_atom_types.size() < 1) {
			sprintf(err_msg, "No atoms found in system");
			flag_error (err_msg, __FILE__, __LINE__);
			return ILLEGAL_VALUE;
		}
		
		// check integrator has been set
		if (integrator == NULL) {
			sprintf(err_msg, "No integrator set on rank %d", rank);
			flag_error (err_msg, __FILE__, __LINE__);
			return ILLEGAL_VALUE;
		}
		
		// check timesteps, dt > 0
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
		
		if (sys->max_rcut() > 0.5*sys->box()[0] / nprocs) {
			sprintf(err_msg, "The maximum r_cut is too large for the system. Use a smaller rcut or fewer processors.");
			flag_error (err_msg, __FILE__, __LINE__);
			return ILLEGAL_VALUE;
		}

		// setup neighbouring domain info on each processor
		//sys->set_rank (rank);
		
		int print_step;
		if (timesteps < 100) {
			print_step = 1;
		} else {
			print_step = (int) floor(timesteps/100.0);
		}

		// execute loops
		for (int i = 0; i < timesteps; ++i) {
			MPI_Barrier(MPI_COMM_WORLD);
			// move atoms that have left the simulation box
			if (nprocs > 1) {
				check = send_atoms(sys);
				if (check != 0) {
					sprintf(err_msg, "Error encountered during sending atoms after step %d", i+1);
					flag_error (err_msg, __FILE__, __LINE__);
					return check;
				}
			}

			// clear out the forces on each atom which is NECESSARY before each new step
			for (int j = 0; j < sys->natoms(); ++j) {
				for (int k = 0; k < NDIM; ++k) {
					sys->get_atom(j)->force[k] = 0.0;
				}
			}
			
			// calc_force
			check = force_calc(sys);

			if (check != 0) {
				sprintf(err_msg, "Error encountered during force calc after step %d", i+1);
				flag_error (err_msg, __FILE__, __LINE__);
				return check;
			}
			

			// create animation
			write_xyz (outfile, sys, i, wrap_pos);
			
			// step forward
			check = integrator->step(sys);
			if (check != 0) {
				sprintf(err_msg, "Error encountered during integration after step %d", i+1);
				flag_error (err_msg, __FILE__, __LINE__);
				return check;
			}
			
			// report progress
			if (rank == 0) {
				if (i%print_step == 0) {
					sprintf(err_msg, "Finished %d of %d total steps", i, timesteps);
					flag_notify (err_msg, __FILE__, __LINE__);
				}
			}
		
		}
		write_xyz (outfile, sys, timesteps, wrap_pos);
		
	return 0;
	}
}
