/**
 \brief MD Integrator(s) Information
 \author Nathan A. Mahynski
 \file integrator.h
 **/

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "system.h"
#include "misc.h"
#include "mpi.h"
#include "domain_decomp.h"
#include "atom.h"
#include "global.h"
#include "force_calc.h"
#include "read_xml.h"

using namespace sim_system; 
using namespace misc;
using namespace atom;
using namespace std;

//! \namespace integrator
//! Namespace for the base class for integrators
namespace integrator {
	//!< Abstract base class for integrators
	/*!
	 The step function should return 0 if successfully executed, integer error flag otherwise.
	 */
	class Integrator {
	public:
		Integrator(){};
		~Integrator(){};
		void set_dt (const double dt) {dt_ = dt;}
		double dt () const {return dt_;}
		void set_temp (double temp) {temp_ = temp;}
		double getTemp() const {return temp_;}
		virtual int step (System *sys) = 0;		//!< Requires all subclasses to be able to execute a step
		double type();
	protected:
		double dt_;
		double temp_;
	};
	
	//! NVE, Verlet
	class Verlet : public Integrator {
	public:
		Verlet(double deltat);
		~Verlet(){};
		double getTime();
		double getdt();
		int step(System *sys);
		double type; 
	private:
		int timestep_;
		double dt2_;
	};
	
	//! Loop through all atoms in the system, send atoms that are outside this domain to the proper domain and delete them locally; also receive and insert new ones.
	int move_atoms (System *sys, const int rank, const int nprocs);
	
	//! NVE, Velocity Verlet
	class Velocity_verlet : public Integrator {
	public:
	private:
	};
  
	//! NVT, Andersen thermostat
	class Andersen : public Integrator {
	public:
		Andersen(double deltat, double temp);
		~Andersen(){};
		void set_temp (double temp) {temp_ = temp;}
		double getTemp() const {return temp_;}
		double getdt();
		int step(System *sys);
		double type;
	private:
		int timestep_;
		double dt2_;
		double temp_;
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
	int run (System *sys, Integrator *integrator, const int timesteps, const string outfile);
}


#endif
