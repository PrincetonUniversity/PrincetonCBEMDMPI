/**
 \file read_xml.h
 \brief I/O for XML
 \author Nathan A. Mahynski
 **/

#ifndef READ_XML_H_
#define READ_XML_H_

#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include "mpi.h"
#include "misc.h"
#include "atom.h"
#include "system.h"
#include "global.h"
#include "domain_decomp.h"
#include "read_interaction.h"

using namespace std;
using namespace sim_system;
using namespace misc;
using namespace atom;
using namespace boost::algorithm;

int read_xml (const string filename, System *sys);
int print_xml (const string filename, const System *sys);
int write_xyz (const string filename, const System *sys, const int timestep, const bool wrap_pos);

#endif
