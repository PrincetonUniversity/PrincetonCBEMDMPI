/*!
 \brief Initialization routines for System object.
 \author Nathan A. Mahynski
 */

#ifndef INITIALIZE_H_
#define INITIALIZE_H_

#include "read_xml.h"
#include "read_interaction.h"
#include "global.h"
#include "mpi.h"
#include "system.h"

using namespace std;
using namespace sim_system;
using namespace misc;
using namespace atom;

//! Parse an XML and energy file to obtain atom and interaction information. 
int initialize_from_files (const string xml_filename, const string energy_filename, System *sys);
//! Call MPI_Init and start the MPI ensuring it began successfully.
int start_mpi (int argc, char *argv[]);
//! Finalize MPI
int end_mpi ();
//! Abort MPI
int abort_mpi ();

#endif
