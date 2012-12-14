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

int read_interactions(const string filename, System *sys);

force_energy_ptr get_fn(const string name);

#endif
