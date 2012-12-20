#include "force_calc.h"

using namespace std;
using namespace sim_system;

int force_calc(System *sys) { // pass interaction array
	const vector<double> box = sys->box();
	double kinetic_energy = 0.0, potential_energy = 0.0, dE, totKE, totPE;
	for (int i=0; i!=sys->total_atoms(); ++i) {
		for (int j=i+1; j!=sys->total_atoms(); ++j) {
		    if ( (*(sys->get_atom(i))).sys_index != (*(sys->get_atom(j))).sys_index ) {
				try {
					dE = sys->interact[(*(sys->get_atom(i))).sys_index][(*(sys->get_atom(j))).sys_index].force_energy(sys->get_atom(i), sys->get_atom(j), &box);
				}
				catch (exception& e) {
					flag_error(e.what(), __FILE__, __LINE__);
					return ILLEGAL_VALUE;
				}
				
				// don't count ghosts
				if (i < sys->natoms() && j < sys->natoms()) {
					potential_energy += dE;
					for (int k = 0; k < 3; ++k) {
						kinetic_energy += 0.5*(sys->get_atom(i)->mass*sys->get_atom(i)->vel[k]*sys->get_atom(i)->vel[k]+sys->get_atom(j)->mass*sys->get_atom(j)->vel[k]*sys->get_atom(j)->vel[k]);
					}
				}
		    }
		}
	}
	
	// keep track of these on all processors (needed for things like thermostats, etc.)
	MPI_Allreduce (&kinetic_energy, &totKE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce (&potential_energy, &totPE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	sys->set_total_KE(totKE);
	sys->set_total_PE(totPE);
	return 0;
}
