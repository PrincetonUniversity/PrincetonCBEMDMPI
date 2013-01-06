/**
 \brief MD Integrator(s) Information
 \author {Nathan A. Mahynski, George Khoury}
 \file integrator.h
 **/

#ifndef INTEGRATOR_H_
#define INTEGRATOR_H_

#include "system.h"
#include "misc.h"
#include "mpi.h"
#include "atom.h"
#include "global.h"
#include "force_calc.h"
#include "read_xml.h"
#include "force_calc.h"
#include <boost/tr1/random.hpp>

using namespace std;

//!< Abstract base class for integrators
/*!
 The step function should return SAFE_EXIT if successfully executed, integer error flag otherwise.
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
	
private:		
	int timestep_;
	double dt2_;							// dt^2
};

//! NVT, Andersen thermostat
class Andersen : public Integrator {
public:
	Andersen(double deltat, double temp, double nu);
	~Andersen(){};
	void set_temp (double temp) {temp_ = temp;}
	double getTemp() const {return temp_;}
	double getdt();
	int step(System *sys);
	double nu () const {return nu_;}
	
private:
	int timestep_;							// Current timestep needs to be recorded for Andersen
	double dt2_;
	double temp_;
	double nu_;
};
  	
//!< Run (i.e. integrate) a system forward in time for a specified number of timesteps
int run (System *sys, Integrator *integrator, const int timesteps, const string outfile);

#endif
