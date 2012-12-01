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
	  Verlet(double deltat);
	  ~Verlet();
	  double getTime();
	  double getdt();
	  int step(System *sys);
	private:
	  int timestep_;
	  double dt_;
	  double dt2_;
	  vector <double[3]> prev_pos_;
	};
	
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
	int run (System *sys, const Integrator *integrator, const int timeteps) ;

}


#endif
