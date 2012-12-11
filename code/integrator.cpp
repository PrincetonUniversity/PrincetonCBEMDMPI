/**
   MD Integrator(s) Information
   \authors{George Khoury, Carmeline Dsilva, Nathan A. Mahynski}
**/

#include "integrator.h"

namespace integrator {
	//! NVE, Verlet
	Verlet::Verlet (double deltat) {
		timestep_ = 0;
		dt_ = deltat;
		dt2_ = dt_ * dt_;
	}

	// example:
	/*! 
	  This function integrates.
	  Call it to make a step.
	  Will update atom positions and velocities
	*/
	int Verlet::step (System *sys) {
		// write the code to update positions, etc. here
		// you can plan on being able to call a calc_force routine that calculates and stores the instantaneous force in
		// the cartesian directions
		/* Arun Prabhu: The next four lines allocate the appropriate ammount of memory
		   for the variable prev_pos which is a vector of vectors but must represent
		   a matrix of size m x n where m = number of atoms in system,
		   and n = 3 (number of spatial dimensions;
		   This comment may be removed later once everyone is aware of the issue and the fix */
		prev_pos_.resize(sys->natoms());
		for (int i=0; i<sys->natoms(); i++) {
		  prev_pos_[i].reserve(3);
		}
		double prev_prev_pos;
		if (timestep_ == 0) {
			for (int i = 0; i < sys->natoms(); ++i) {
				for (int j = 0; j < 3; ++j) {
					prev_pos_[i][j] = sys->get_atom(i)->pos[j];
					sys->get_atom(i)->pos[j] += sys->get_atom(i)->vel[j] * dt_ + 0.5 * sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
					sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - prev_pos_[i][j]) / dt_;
				}
			}
		}
		else {
			for (int i = 0; i < sys->natoms(); ++i) {
				for (int j = 0; j < 3; ++j) {
					prev_prev_pos = prev_pos_[i][j];
					prev_pos_[i][j] = sys->get_atom(i)->pos[j];
					sys->get_atom(i)->pos[j] = 2.0 *  prev_pos_[i][j] - prev_prev_pos + sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
					sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - prev_prev_pos) / (2.0 * dt_);
				}
			}
		}
		timestep_++;
		return 0;
	}
	
	//! Loop through all atoms in the system, send atoms that are outside this domain to the proper domain.
	int move_atoms (System *sys, const int rank, const int nprocs) {
		int *nsend_atoms = (int *) calloc (nprocs, sizeof(int));
		int *nrecv_atoms = (int *) calloc (nprocs, sizeof(int));
		int iproc;
		vector <double> pos(3);
		MPI_Request *reqs_nsend = (MPI_Request *) malloc(nprocs*sizeof(MPI_Request)), *reqs_nrecv = (MPI_Request *) malloc(nprocs*sizeof(MPI_Request));
		MPI_Status *stat_recv = (MPI_Status *) malloc(nprocs*sizeof(MPI_Status));
		Atom *leaving_atoms[nprocs], *arriving_atoms[nprocs];
		vector < vector <int> > atom_goes;
		vector <int> atom_goes_loc;
		
		try {
			atom_goes.resize(nprocs);
		}
		catch (bad_alloc& ba) {
			sprintf(err_msg, "Couldn't allocate space for atom send list");
			flag_error (err_msg, __FILE__, __LINE__);
			return -1;
		}
		
		// collect and enumerate number atoms that go to each processor
		int tot_send = 0;
		for (int i = 0; i < sys->natoms(); ++i) {
			for (int j = 0; j < 3; ++j) {
				pos[j] = sys->get_atom(i)->pos[j];
			}
			iproc = get_processor (pos, sys->proc_widths, sys->final_proc_breakup);
			if (iproc != rank) {
				nsend_atoms[iproc]++;
				tot_send++;
				atom_goes[iproc].push_back(i);
				// might be faster to allocate after this look and assign later since no reallocation every time, but this way atoms are ordered...
				atom_goes_loc.push_back(i);
			}
		}

		// send/recv how many are going to/coming from each domain
		for (int i = 0; i < nprocs; ++i) {
			if (i == rank) {
				continue;
			}
			MPI_Isend (&nsend_atoms[i], 1, MPI_INT, i, rank, MPI_COMM_WORLD, &reqs_nsend[i]);
			MPI_Irecv (&nrecv_atoms[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &reqs_nrecv[i]);
		}
		MPI_Waitall (nprocs, reqs_nrecv, stat_recv);
		
		// allocate atoms
		for (int i = 0; i < nprocs; ++i) {
			if (i == rank) {
				continue;
			}
			leaving_atoms[i] = new Atom[nsend_atoms[i]];
			arriving_atoms[i] = new Atom[nrecv_atoms[i]];
		}
		
		// assign leaving atoms ??
		for (int i = 0; i < nprocs; ++i) {
			if (i != rank) {
				for (int j = 0; j < nsend_atoms[i]; ++j) {
					// create copy to send to new processor
					leaving_atoms[i][j] = sys->copy_atom(atom_goes[i][j]);
				}
			}
		}
		
		// delete originals from here
		sys->delete_atoms(atom_goes_loc);
		
		// send/receive atoms themselves from all domains
		for (int i = 0; i < nprocs; ++i) {
			if (i == rank) {
				continue;
			} 
			if (nrecv_atoms[i] > 0) {
				MPI_Irecv (&arriving_atoms[i], nrecv_atoms[i], MPI_ATOM, i, i, MPI_COMM_WORLD, &reqs_nrecv[i]);
			} else {
				reqs_nrecv[i] = MPI_SUCCESS;
			}

			if (nsend_atoms[i] > 0) {
				MPI_Isend (&leaving_atoms[i], nsend_atoms[i], MPI_ATOM, i, rank, MPI_COMM_WORLD, &reqs_nsend[i]);
			} else {
				// this is probably superfluous
				reqs_nsend[i] = MPI_SUCCESS;
			}
		}
		MPI_Waitall (nprocs, reqs_nrecv, stat_recv);
		
		// add ones that arrived
		for (int i = 0; i < nprocs; ++i) {
			if (i == rank || nrecv_atoms[i] == 0) {
				continue;
			} else {
				sys->add_atoms(nrecv_atoms[i], arriving_atoms[i]);
			}
		}
		
		for (int i = 0; i < nprocs; ++i) {
			if (i == rank) {
				continue;
			}
			delete [] leaving_atoms[i]; 
			delete [] arriving_atoms[i];
		}
		
		// other deallocs
		free (reqs_nsend);
		free (reqs_nrecv);
		free (stat_recv);
		free (nsend_atoms);
		free (nrecv_atoms);
		
		return SAFE_EXIT;
	}
	
	//!< Run (i.e. integrate) a system forward in time for a specified number of timesteps
	/*!
	  Function returns 0 if successful, -1 if it encountered an error.
	  \param [in] \*sys Pointer to System to integrate
	  \param [in] \*sys Pointer to Integrator to use
	  \param [in] timesteps Number of timesteps to integrate over
	*/
	int run (System *sys, Integrator *integrator, const int timesteps) {
		// The way this function is written it can be easily interpreted by SWIG with python!
		int check = 0, nprocs, rank;
		char err_msg[MYERR_FLAG_SIZE]; //!< Error message buffer commonly used in routines in this namespace
  
		MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
		MPI_Comm_rank (MPI_COMM_WORLD, &rank);
		
		// before starting, need to check that all requisite variables are set
  
		
		// execute loops
		MPI_Barrier(MPI_COMM_WORLD);
		for (int i = 0; i < timesteps; ++i) {
			// calc_force
			
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			// step forward
			check = integrator->step(sys);
			if (check != 0) {
				sprintf(err_msg, "Error encountered during integration after step %d", i+1);
				flag_error (err_msg, __FILE__, __LINE__);
				return check;
			}
			
			// move atoms as necessary
			move_atoms (sys, rank, nprocs);
			
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}


}
