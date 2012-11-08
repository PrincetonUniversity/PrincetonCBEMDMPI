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
		int set_mpi (const int flag);			//!< Set flag to indicate if using MPI
		int set_gpu (const int flag);			//!< Set flag to indicate if using GPU's
		
	private:
		int use_mpi_;							//!< Flag to indicate if using MPI to perform calculations
		int use_gpu_;							//!< Flag to indicate if using GPU's to perform calculations
	};
	
	//! Default constructor sets parallel implementation flags to 0 (false) by default
	Integrator::Integrator () {
		use_gpu_ = 0;
		use_mpi_ = 0;
	}
	
	/*!
	 \param [in] flag 0 for False, 1 for True.  Anything else will return an error and leave the flag alone.
	 */
	int Integrator::set_mpi (const int flag) {
		if (flag != 0 && flag != 1) {
			sprintf(err_msg, "Illegal MPI flag for integrator", __FILE__, __LINE__);
			return ILLEGAL_VALUE;
		} else {
			use_mpi_ = flag;
		}
	}
	
	/*!
	 \param [in] flag 0 for False, 1 for True.  Anything else will return an error and leave the flag alone.
	 */
	int Integrator::set_gpu (const int flag) {
		if (flag != 0 && flag != 1) {
			sprintf(err_msg, "Illegal GPU flag for integrator", __FILE__, __LINE__);
			return ILLEGAL_VALUE;
		} else {
			use_gpu_ = flag;
		}
	}
	
	// More specific integrators include:
	
	//! NVE, Verlet
	class Verlet : public Integrator {
	public:
	private:
	};
	
	// example:
	int Verlet::step (System *sys) {
		// write the code to update positions, etc. here
		// (but this should really go in a .cc file)
		
	}
	
	//! NVE, Velocity Verlet
	class Velocity_verlet : public Integrator {
	public:
	private:
	};
	
	//! NVE, Beeman algorithm
	class Beeman : public Integrator {
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
		// Inside this function we should initially check if the system is "prepared", i.e. all necessary vars specified
		
		// loop for timestep steps, report any errors as encountered
		
		// That's it!
	}
}


#endif
