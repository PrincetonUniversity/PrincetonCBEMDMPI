/**
   \brief MD Integrator(s) Information
   \authors{George Khoury, Carmeline Dsilva, Nathan A. Mahynski}
   \file integrator.cpp
**/

#include "integrator.h"
#include "force_calc.h"
#include </home/gkhoury/boost_1_52_0/boost/tr1/random.hpp>


//! \namespace integrator
//! Namespace for the base class for integrators
namespace integrator {

	//
	// Generate a random number between 0 and 1
	// return a uniform number in [0,1].
	double unifRand()
	{
	    return rand() / double(RAND_MAX);
	}
	//
	// Generate a random number in a real interval.
	// param a one end point of the interval
	// param b the other end of the interval
	// return a inform rand numberin [a,b].
	double unifRand(double a, double b)
	{
	    return (b-a)*unifRand() + a;
	}


	// ! NVT, Andersen Thermostat
	Andersen::Andersen(double deltat, double temp){
		timestep_ = 0;
		dt_ = deltat;
		dt2_ = dt_*dt_;
		temp_ = temp;
		type = 2;
	}
	
	// http://phycomp.technion.ac.il/~talimu/md4.html
	int Andersen::step (System *sys) {
		double prev_prev_pos;
		double nu = 100.0; // thermostat parameter controlling frequency of collisions with heat bath
		vector <double> box = sys->box();
		int check = 0, nprocs, rank;
		char err_msg[MYERR_FLAG_SIZE];
		MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
		MPI_Comm_rank (MPI_COMM_WORLD, &rank);
                //std::tr1::default_random_engine generator;
		double totalmass = 0;  // initialize totalmass of system
		for (int i = 0; i < sys->natoms(); ++i) {
			totalmass += sys->get_atom(i)->mass;
		}
		//totalmass = totalmass / sys->natoms();
		//std::mt19937 Myceng;
		typedef std::tr1::linear_congruential<int, 16807, 0, (int)((1U << 31) -1 ) > Myceng;
		Myceng eng;
                cout << "temp is " << temp_ << endl;
//		if (timestep_ == 0) {
			// step #1
			for (int i = 0; i < sys->natoms(); ++i) {
				for (int j = 0; j < NDIM; ++j) {
					sys->get_atom(i)->prev_pos[j] = sys->get_atom(i)->pos[j];
					sys->get_atom(i)->pos[j] += sys->get_atom(i)->vel[j] * dt_ + 0.5 * sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
					sys->get_atom(i)->vel[j] += 0.5* dt_*sys->get_atom(i)->force[j] / sys->get_atom(i)->mass;
					//sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - sys->get_atom(i)->prev_pos[j]) / dt_;
				}
			}

			// force calc
			check = force_calc(sys);

			// move atoms
			// check to move atoms if necessary
			if (nprocs > 1) {
				check = move_atoms (sys, rank, nprocs);
				if (check != 0) {
					sprintf(err_msg, "Couldn't move atoms from rank %d on step %d", rank, dt_);
					flag_error (err_msg, __FILE__, __LINE__);
					return check;
				}
			}

			// step #2 
			double tempa = 0;
			//temp_ = 0 ;
			for (int i = 0; i < sys->natoms(); ++i) {
				for (int j = 0; j < NDIM; ++j) {
					prev_prev_pos = sys->get_atom(i)->prev_pos[j];
					sys->get_atom(i)->prev_pos[j] = sys->get_atom(i)->pos[j];
					//sys->get_atom(i)->pos[j] = 2.0 *  sys->get_atom(i)->prev_pos[j] - prev_prev_pos + sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
					sys->get_atom(i)->vel[j] = sys->get_atom(i)->vel[j] + 0.5 * dt_ * sys->get_atom(i)->force[j] / sys->get_atom(i)->mass;
					tempa += sys->get_atom(i)->mass*sys->get_atom(i)->vel[j]*sys->get_atom(i)->vel[j];
					//temp_ += sys->get_atom(i)->mass*sys->get_atom(i)->vel[j]*sys->get_atom(i)->vel[j];
				}
			}

			double s = 3.0 ;
			tempa = tempa / ( s * (sys->natoms()));
			//temp_ = temp_ / (s* sys->natoms());

			cout << "Instantaneous Temp is "  << tempa << endl;
			//cout << "Instantaneous Temp is " << temp_ << endl;
			double sig;
			sig = sqrt(temp_); //should this be div by totalmass?!
			std::tr1::normal_distribution<double> distribution(0.0,sig);
			double rannum;
			for (int i = 0; i < sys->natoms(); ++i) {
				rannum = unifRand();
				//cout << rannum << " is random number " << endl;
				if (rannum < nu*dt_) {
					//sig = sqrt(temp_/sys->get_atom(i)->mass);
					//std::tr1::normal_distribution<double> distribution(0.0,sig);
					for (int j = 0; j < NDIM; ++j) {
					//if (unifRand() < nu*dt_) {
						sys->get_atom(i)->vel[j] = distribution(eng);
					}
				}
			}
//		}
//		else {

/*			// step 1
			for (int i = 0; i < sys->natoms(); ++i) {
				for (int j = 0; j < NDIM; ++j) {
					prev_prev_pos = sys->get_atom(i)->prev_pos[j];
					sys->get_atom(i)->prev_pos[j] = sys->get_atom(i)->pos[j];
					sys->get_atom(i)->pos[j] = 2.0 *  sys->get_atom(i)->prev_pos[j] - prev_prev_pos + sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
                                        sys->get_atom(i)->vel[j] += 0.5* dt_*sys->get_atom(i)->force[j] / sys->get_atom(i)->mass;

					//sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - sys->get_atom(i)->prev_pos[j]) / dt_;
				}
			}

                        // force calc
                        check = force_calc(sys);

			// move atoms
			// check to move atoms if necessary
			if (nprocs > 1) {
				check = move_atoms (sys, rank, nprocs);
				if (check != 0) {
					sprintf(err_msg, "Couldn't move atoms from rank %d on step %d", rank, dt_);
					flag_error (err_msg, __FILE__, __LINE__);
					return check;
				}
			}

			// step 2
			double tempa = 0;
			//temp_ = 0;
			for (int i = 0; i < sys->natoms(); ++i) {
				for (int j = 0; j < NDIM; ++j) {
					prev_prev_pos = sys->get_atom(i)->prev_pos[j];
					sys->get_atom(i)->prev_pos[j] = sys->get_atom(i)->pos[j];
					sys->get_atom(i)->pos[j] = 2.0 *  sys->get_atom(i)->prev_pos[j] - prev_prev_pos + sys->get_atom(i)->force[j] / sys->get_atom(i)->mass * dt2_;
					//sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - sys->get_atom(i)->prev_pos[j]) / dt_;
                                        sys->get_atom(i)->vel[j] = sys->get_atom(i)->vel[j] + 0.5 * dt_ * sys->get_atom(i)->force[j];
                                        tempa += sys->get_atom(i)->mass*sys->get_atom(i)->vel[j]*sys->get_atom(i)->vel[j];
					//temp_ += sys->get_atom(i)->mass*sys->get_atom(i)->vel[j]*sys->get_atom(i)->vel[j];
				}
			}

                        double s = 3.0 / totalmass;
                        tempa = tempa / ( (s * sys->natoms()));
			//temp_ = temp_ / (s* sys->natoms());

			cout << "Instantaneous Temp " << tempa << endl;
			//cout << "Instantaneous Temp " << temp_ << endl;
			double sig;
			sig = sqrt(temp_); // should this be tempa NO
			std::tr1::normal_distribution<double> distribution(0.0,sig);
			for (int i = 0; i < sys->natoms(); ++i) {
				for (int j = 0; j < NDIM; ++j) {
					if (unifRand() < nu*dt_) {
						sys->get_atom(i)->vel[j] = distribution(eng);
					}
				}
			}
//		}
		*/
		timestep_++;
		return 0;
	}
	//! NVE, Verlet
	Verlet::Verlet (double deltat) {
		timestep_ = 0;
		dt_ = deltat;
		dt2_ = dt_ * dt_;
		type = 1;
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
						}
					*/
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
					sys->get_atom(i)->vel[j] = (sys->get_atom(i)->pos[j] - sys->get_atom(i)->prev_pos[j]) / dt_;
				}
			}
		}
		timestep_++;
		return 0;
	}
	
	
	/*! 
	  This function loops through all atoms in the system, sends atoms that
      are outside this domain to the proper domain and
      deletes them locally. Also, it will receive and insert new atoms.
   	  \param [in] sys Pointer to System to integrate
	  \param [in] rank processor rank
	  \param [in] nprocs number of processors 
	*/	
	int move_atoms (System *sys, const int rank, const int nprocs) {
	    char err_msg[MYERR_FLAG_SIZE];
		vector <int> nsend_atoms(nprocs,0), nrecv_atoms(nprocs,0);
		int iproc;
		vector <double> pos(NDIM);
		vector < vector <int> > atom_goes;
		vector <int> atom_goes_loc, dummy_catch;
		Atom **leaving_atoms = (Atom **)malloc(nprocs*sizeof(Atom *)), **arriving_atoms = (Atom **)malloc(nprocs*sizeof(Atom *));
		const int num_procs=nprocs;
		MPI_Request reqs_nsend[num_procs], reqs_nrecv[num_procs], reqs_nrecv2[num_procs];
		MPI_Status stat_recv[num_procs];

		try {
			atom_goes.resize(nprocs);
		}
		catch (bad_alloc& ba) {
			sprintf(err_msg, "Couldn't allocate space for atom send list");
			flag_error (err_msg, __FILE__, __LINE__);
			return BAD_MEM;
		}

		// collect and enumerate number atoms that go to each processor
		int tot_send = 0;
		for (int i = 0; i < sys->natoms(); ++i) {
			for (int j = 0; j < NDIM; ++j) {
				pos[j] = sys->get_atom(i)->pos[j];
			}
			iproc = get_processor (pos, sys);
			if (iproc < 0 || iproc >= nprocs) {
				sprintf(err_msg, "Processor %d out of bounds for atom %d, coordinates = (%g, %g, %g)", iproc, sys->get_atom(i)->sys_index, pos[0], pos[1], pos[2]);
				flag_error (err_msg, __FILE__, __LINE__);
				return ILLEGAL_VALUE;
			}
			if (iproc != rank) {
				nsend_atoms[iproc]++;
				tot_send++;
				atom_goes[iproc].push_back(i);
				// might be faster to allocate after this look and assign later since no reallocation every time, but this way atoms are ordered...
				atom_goes_loc.push_back(i);
			} else {
				// don't bother sending/recv atoms from yourself
				nsend_atoms[iproc] = 0;
			}
		}

		// send/recv how many are going to/coming from each domain (including self technically, but overhead is very low)
		for (int i = 0; i < nprocs; ++i) {
			MPI_Isend (&nsend_atoms[i], 1, MPI_INT, i, rank, MPI_COMM_WORLD, &reqs_nsend[i]);
			MPI_Irecv (&nrecv_atoms[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, &reqs_nrecv[i]);
		}
		MPI_Waitall (nprocs, reqs_nrecv, stat_recv);

		// allocate atoms
		for (int i = 0; i < nprocs; ++i) {
			if (nsend_atoms[i] > 0) {
				leaving_atoms[i] = (Atom *)malloc(nsend_atoms[i]*sizeof(Atom));
			} else {
				leaving_atoms[i] = (Atom *)malloc(1*sizeof(Atom));
			}
			if (nrecv_atoms[i] > 0) {
				arriving_atoms[i] = (Atom *)malloc(nrecv_atoms[i]*sizeof(Atom));
			} else {
				arriving_atoms[i] = (Atom *)malloc(1*sizeof(Atom));
			}
		}

		// assign leaving atoms
		for (int i = 0; i < nprocs; ++i) {
			if (i != rank) {
				for (int j = 0; j < nsend_atoms[i]; ++j) {
					leaving_atoms[i][j] = sys->copy_atom(atom_goes[i][j]);
				}
			}
		}

		// delete originals from the local system
		sys->delete_atoms(atom_goes_loc);
		
		// send/receive atoms themselves from all domains
		for (int i = 0; i < nprocs; ++i) {
			MPI_Irecv (&arriving_atoms[i][0], nrecv_atoms[i], MPI_ATOM, i, i, MPI_COMM_WORLD, &reqs_nrecv2[i]);
			MPI_Isend (&leaving_atoms[i][0], nsend_atoms[i], MPI_ATOM, i, rank, MPI_COMM_WORLD, &reqs_nsend[i]);
		}
		MPI_Waitall (nprocs, reqs_nrecv2, stat_recv);

		// add ones that arrived
		for (int i = 0; i < nprocs; ++i) {
			if (i == rank || nrecv_atoms[i] == 0) {
				continue;
			} else {
				dummy_catch = sys->add_atoms(nrecv_atoms[i], &arriving_atoms[i][0]);
			}
		}

		for (int i = 0; i < nprocs; ++i) {
			free (leaving_atoms[i]);
			free (arriving_atoms[i]);
		}
		free (leaving_atoms);
		free (arriving_atoms);
		return SAFE_EXIT;
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
		
		// setup neighbouring domain info on each processor
		sys->set_rank (rank);
		if (sys->gen_domain_info()) {
		    sprintf(err_msg, "Problem generating basic info for processor rank %d", rank);
		    flag_error (err_msg, __FILE__, __LINE__);
		    return ILLEGAL_VALUE; //This needs to be corrected, illegal value is here temporarily
		}
		if (gen_send_table(sys)) {
		    sprintf(err_msg, "Problem generating table of neighbouring processors on rank %d", rank);
		    flag_error (err_msg, __FILE__, __LINE__);
		    return ILLEGAL_VALUE; //This needs to be corrected, illegal value is here temporarily
		}
		
		int print_step;
		if (timesteps < 100) {
			print_step = 1;
		} else {
			print_step = (int) floor(timesteps/100.0);
		}

		// execute loops
		MPI_Barrier(MPI_COMM_WORLD);
		for (int i = 0; i < timesteps; ++i) {
			// clear out the forces on each atom which is NECESSARY before each new step
			for (int j = 0; j < sys->natoms(); ++j) {
				for (int k = 0; k < NDIM; ++k) {
					sys->get_atom(j)->force[k] = 0.0;
				}
			}
			
		    // generate lists of atoms to be sent to neighbouring cells
		    if (gen_send_lists(sys)) {
				sprintf(err_msg, "Problem generating lists to send out on rank %d", rank);
				flag_error (err_msg, __FILE__, __LINE__);
				return ILLEGAL_VALUE; //This needs to be corrected, illegal value is here temporarily
		    }
			
		    // communicate the above atoms (ghost atoms) to the appropriate processor
		    if (communicate_skin_atoms(sys)) {
				sprintf(err_msg, "Problem communicating lists on rank %d", rank);
				flag_error (err_msg, __FILE__, __LINE__);
				return ILLEGAL_VALUE; //This needs to be corrected, illegal value is here temporarily
		    }
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			// calc_force
			if (i == 0) {
				check = force_calc(sys);
			}
			//else if (i > 0 and integrator->type == 1){
			else if (integrator->getTemp() == 0 ) {
				cout << " You are running NVE,nothing changed " << endl;
				check = force_calc(sys);
			}
			else {
				cout << "You are running NVT, nothing changed " << endl;
			}

			if (check != 0) {
				sprintf(err_msg, "Error encountered during force calc after step %d", i+1);
				flag_error (err_msg, __FILE__, __LINE__);
				return check;
			}
			
			// Delete the atoms the processor is not responsible for
			sys->clear_ghost_atoms ();


			// create animation
			write_xyz (outfile, sys, i, wrap_pos);

			// step forward
			check = integrator->step(sys);
			if (check != 0) {
				sprintf(err_msg, "Error encountered during integration after step %d", i+1);
				flag_error (err_msg, __FILE__, __LINE__);
				return check;
			}
			
			// check to move atoms if necessary
			if (nprocs > 1) {
				check = move_atoms (sys, rank, nprocs);
				if (check != 0) {
					sprintf(err_msg, "Couldn't move atoms from rank %d on step %d", rank, i+1);
					flag_error (err_msg, __FILE__, __LINE__);
					return check;
				}
			}
			
			MPI_Barrier(MPI_COMM_WORLD);
			
			// report progress
			if (rank == 0) {
				if (i%print_step == 0) {
					sprintf(err_msg, "Finished %d of %d total steps", i, timesteps);
					flag_notify (err_msg, __FILE__, __LINE__);
				}
			}

		}
		write_xyz ("output.xyz", sys, timesteps, wrap_pos);
		
	return 0;
	}
}
