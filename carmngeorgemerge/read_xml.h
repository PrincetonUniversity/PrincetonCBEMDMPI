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
#include "read_interaction.h"

using namespace std;
using namespace boost::algorithm;

//! Read an xml file to initialize a System object
int read_xml (const string filename, System *sys);

//! Print the current state of the system to a .xml file
int print_xml (const string filename, const System *sys);

//! Write an animation file (.xyz)
int write_xyz (const string filename, const System *sys, const int timestep, const bool wrap_pos);

#endif
