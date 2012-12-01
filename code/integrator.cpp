/**
   MD Integrator(s) Information
   \authors{George Khoury, Carmeline Dsilva}
**/

#include "integrator.h"

  //! NVE, Verlet
  Verlet::Verlet(double deltat) {
    timestep = 0;
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
    prev_pos_.resize(sys.natoms());
    if (time == 0) {
      for (int i = 0; i < sys.natoms(); ++i) {
	for (int j = 0; j < 3; ++j) {
	  prev_pos_[i][j] = sys.get_atom(i).pos[j];
	  sys.get_atom(i).pos[j] += sys.get_atom(i).vel[j] * dt_ + 0.5 * sys.get_atom(i).force[j] / sys.get_atom(i).mass * dt2_;
	  sys.get_atom(i).vel[j] = (sys.get_atom(i).pos[j] - prev_pos_[i][j]) / dt_;
	}
      }
    }
  
    double prev_prev_pos;
    else {
      for (int i = 0; i < sys.natoms(); ++i) {
	for (int j = 0; j < 3; ++j) {
	  prev_prev_pos = prev_pos_[i][j];
	  prev_pos_[i][j] = sys.get_atom(i).pos[j];
	  sys.get_atom(i).pos[j] = 2.0 *  prev_pos_[i][j] - prev_prev_pos + sys.get_atom(i).force[j] / sys.get_atom(i).mass * dt2_;
	  sys.get_atom(i).vel[j] = (sys.get_atom(i).pos[j] - prev_prev_pos) / (2.0 * dt_);
	}
      }
    }
    timestep++;
  }
	
  //!< Run (i.e. integrate) a system forward in time for a specified number of timesteps
  /*!
    Function returns 0 if successful, -1 if it encountered an error.
    \param [in] \*sys Pointer to System to integrate
    \param [in] \*sys Pointer to Integrator to use
    \param [in] timesteps Number of timesteps to integrate over
  */
  int run (System *sys, const Integrator *integrator, const int timeteps) {
    // The way this function is written it can be easily interpreted by SWIG with python!
    int check = 0;
    
    // before starting, need to check that all requisite variables are set
    
    // execute loops
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = 0; i < timeteps; ++i) {
      check = integrator->step(sys);
      if (check != 0) {
	sprintf(err_msg, "Error encountered during integration after step %d", i+1);
	flag_error (err_msg, __FILE__, __LINE__);
	return check;
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }


