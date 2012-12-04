/**
 I/O for XML
 \author Nathan A. Mahynski
 **/

#ifndef READ_XML_H_
#define READ_XML_H_

#include "system.h"
#include <boost/algorithm/string.hpp>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include "mpi.h"
#include "misc.h"
#include "atom.h"

using namespace std;
using namespace sim_system;
using namespace misc;
using namespace atom;

int read_xml (const string filename, const int nprocs, System *sys);
int print_xml (const string filename, const int nprocs, const int rank, MPI_COMM *MCOMM, System *sys)

#endif
