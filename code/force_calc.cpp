/*!
 \file force_calc.cpp
 \brief Source code for force calculation function
 \authors{Nathan Mahynski, Carmeline Dsilva, Arun Prabhu, George Khoury}
**/

#include "force_calc.h"

using namespace std;
using namespace sim_system;

/*!
 \brief Calculates the forces between the particles in the system
 \param sys [in] System for which to evaluate the forces
*/
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
				
				// don't count ghost-ghost interactions
				// count ghost-cell interactions if the index of the atom in the cell is smaller than the ghost index
				if ((i < sys->natoms() && j < sys->natoms()) || 
				    (i < sys->natoms() && sys->get_atom(i)->sys_index < sys->get_atom(j)->sys_index) ||
				    (j < sys->natoms() && sys->get_atom(j)->sys_index < sys->get_atom(i)->sys_index)) {
					potential_energy += dE;
				}
		    }
		}
		// KE = sum(i,1/2 *m(i)*v(i)*v(i))
		if (i < sys->natoms()) {
			for (int k = 0; k < 3; ++k) {
				kinetic_energy += 0.5*(sys->get_atom(i)->mass*sys->get_atom(i)->vel[k]*sys->get_atom(i)->vel[k]);
			}
		}
	}
	
	
	// keep track of these on all processors (needed for things like thermostats, etc.)
	MPI_Allreduce (&kinetic_energy, &totKE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce (&potential_energy, &totPE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	sys->set_total_KE(totKE);
	sys->set_total_PE(totPE);
	if (sys->rank() == 0) {
		double totE = totPE + totKE; 
		cout << "KE = " << totKE<< ", PE = " << totPE << ", total = " << totE <<  endl;
	}
	return 0;
}
