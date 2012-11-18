/**
 MD Integrator(s) Information
 \author Nathan A. Mahynski
 **/

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "misc.h"
#include "atom.h"
#include "system.h"

using namespace atom;
using namespace std;
using namespace misc;
using namespace system;

namespace integrator {
	char err_msg[ERR_FLAG_SIZE]; //!< Error message buffer commonly used in routines in this namespace
	
	//!< Abstract base class for integrators
	/*!
	 The step function should return 0 if successfully executed, integer error flag otherwise.
	 */
	class Integrator {
	public:
		Integrator();
		~Integrator();
		virtual int step (System *sys) = 0;		//!< Requires all subclasses to be able to execute a step
	};
	
	//! Default constructor sets parallel implementation flags to 0 (false) by default
	Integrator::Integrator () {}
	
	//! NVE, Verlet
	class Verlet : public Integrator {
	public:
	private:
	};
	
	// example:
	int Verlet::step (System *sys) {
		// write the code to update positions, etc. here
		// you can plan on being able to call a calc_force routine that calculates and stores the instantaneous force in
		// the cartesian directions
	}
	
	//! NVE, Velocity Verlet
	class Velocity_verlet : public Integrator {
	public:
	private:
	};
	
	//! NVT, Andersen thermostat
	class Andersen : public Integrator {
	public:
	private:
	};
	
	//! NVT, Nose-Hooser thermostat
	class NoseHoover : public Integrator {
	public:
	private:
	};
	
	//! NPT, Nose-Hoover thermostat + Andersen barostat
	class NHA : public Integrator {
	public:
	private:
	};
	
	
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
}


#endif
