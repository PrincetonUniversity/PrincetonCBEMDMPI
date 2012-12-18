/*!
 Initialization routines for System object.
 \author Nathan A. Mahynski
 */

#include "read_xml.h"
#include "read_interaction.h"
#include "global.h"
#include "mpi.h"
#include "system.h"

using namespace std;
using namespace sim_system;
using namespace misc;
using namespace atom;

int initialize_from_files (const string xml_filename, const string energy_filename, System *sys);
int start_mpi (int argc, char *argv[]);
int end_mpi ();