/*! 
 \file read_interaction.h
 \brief Header file for reading interactions in
**/

#ifndef READ_INTERACTION_H_
#define READ_INTERACTION_H_

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include "mpi.h"
#include "misc.h"
#include "atom.h"
#include "system.h"
#include "global.h"
#include "domain_decomp.h"

using namespace std;
using namespace sim_system;
using namespace misc;
using namespace atom;
using namespace boost::algorithm;

//! Function to read in interaction parameters and store them into the system
int read_interactions(const string filename, System *sys);

//! Function to return a force_energy_ptr given a type of interaction name
force_energy_ptr get_fn(const string name, vector <double> *args, double *r_cut_max);

#endif
