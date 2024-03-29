/*!
 \file force_calc.cpp
 \brief Source code for force calculation
 \authors{Frank Ricci, Jun Park, Nathan Mahynski, Carmeline Dsilva, Arun Prabhu, George Khoury}
**/

#include "force_calc.h"

using namespace std;

/*!
 This function sends the atoms that have left the domain of the processor to the relevant neighboring processor.
 If the atom has moved more than one domain width, it returns an error flag, else returns SAFE_EXIT for success.
 \param \*sys [in] Pointer to system for which to move the atoms from
**/
int send_atoms(System *sys) {
	char err_msg[MYERR_FLAG_SIZE];
	const vector<double> box = sys->box();

	int nprocs, rank;
	MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	
	// store number of atoms to send to and receive from the processor on the left and on the right
	int num_to_left = 0;
	int num_to_right = 0;
	int num_from_left, num_from_right;
	
	// store atoms to send and receive
	vector<Atom> to_left;
	vector<Atom> to_right;
	vector<Atom> from_left;
	vector<Atom> from_right;
	
	// store indices of atoms that have been sent (so we can delete them)
	vector<int> to_delete;

	MPI_Request req[4], req2[4];
	MPI_Status stat[4], stat2[4];

	int proc_to;
	for (int i=0; i!=sys->natoms(); ++i) {
		// calculate the processor for each atom
		proc_to = floor(pbc(sys->get_atom(i)->pos, box)[PARALLELDIM] / box[PARALLELDIM] * nprocs);
		if (proc_to == (rank - 1 + nprocs) % nprocs) {
			to_left.push_back(*(sys->get_atom(i)));
			num_to_left++;
			to_delete.push_back(i);
		}
		else if (proc_to == (rank + 1) % nprocs) {
			to_right.push_back(*(sys->get_atom(i)));
			num_to_right++;
			to_delete.push_back(i);
		}
		else if (proc_to != rank) {
			sprintf(err_msg, "Atom moved too many boxes");
			flag_error (err_msg, __FILE__, __LINE__);
			return ILLEGAL_VALUE;
		}			
	}

	// send number of atoms
	MPI_Isend(&num_to_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req);
	MPI_Irecv(&num_from_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req+1);
	MPI_Isend(&num_to_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+2);
	MPI_Irecv(&num_from_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+3);
	MPI_Waitall (4, req, stat);

	// resize atom vectors
	from_left.resize(num_from_left);
	from_right.resize(num_from_right);
	
	// send atoms
	MPI_Isend(&to_left[0], num_to_left, MPI_ATOM, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2);
	MPI_Irecv(&from_left[0], num_from_left, MPI_ATOM, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2+1);
	MPI_Isend(&to_right[0], num_to_right, MPI_ATOM, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+2);
	MPI_Irecv(&from_right[0], num_from_right, MPI_ATOM, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+3);
	MPI_Waitall (4, req2, stat2);

	// add atoms to system
	sys->add_atoms(&from_left);
	sys->add_atoms(&from_right);

	// delete atoms that we sent to another system
	sys->delete_atoms(to_delete);
	
	MPI_Barrier(MPI_COMM_WORLD);
	return SAFE_EXIT;
}
	
/*!
 Returns SAFE_EXIT if successful, else returns an error flag.
 \param [in] \*sys Pointer to system for which to evaluate the forces
*/
int force_calc(System *sys) { 
	const double skin_cutoff = sys->max_rcut();
	const vector<double> box = sys->box();
	double kinetic_energy = 0.0, potential_energy = 0.0, dE, totKE, totPE;
	int nprocs, rank;
	MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	
	// Store number of atoms to send to and receive from the processor on the left and on the right
	int num_to_left = 0;
	int num_to_right = 0;
	int num_from_left, num_from_right;
	
	// Store atoms to send and receive
	vector<Atom> to_left;
	vector<Atom> to_right;
	vector<Atom> from_left;
	vector<Atom> from_right;

	MPI_Request req[4], req2[4];
	MPI_Status stat[4], stat2[4];

	// Calculate forces between all atoms in the system
	for (int i=0; i!=sys->natoms(); ++i) {
		for (int j=i+1; j!=sys->natoms(); ++j) {
			// A "try" statement is necessary here because interactions can throw errors which need to be caught here
			try {
				dE = sys->interact[(*(sys->get_atom(i))).sys_index][(*(sys->get_atom(j))).sys_index].force_energy(sys->get_atom(i), sys->get_atom(j), &box);
			}
			catch (exception& e) {
				flag_error(e.what(), __FILE__, __LINE__);
				return ILLEGAL_VALUE;
			}
			potential_energy += dE;
		}
		
		// Remember if atom needs to be sent to neighboring processor
		// special case for nprocs==2; don't want to send same atom twise
		if (nprocs == 2) {
 			if (pbc(sys->get_atom(i)->pos, box)[PARALLELDIM] < rank * box[PARALLELDIM] / nprocs + skin_cutoff) {
				to_left.push_back(*(sys->get_atom(i)));
				num_to_left++;
			}
			else if  (pbc(sys->get_atom(i)->pos, box)[PARALLELDIM] > (rank + 1) * box[PARALLELDIM] / nprocs - skin_cutoff) {
				to_right.push_back(*(sys->get_atom(i)));
				num_to_right++;
			}
		} else {
 			if (pbc(sys->get_atom(i)->pos, box)[PARALLELDIM] < rank * box[PARALLELDIM] / nprocs + skin_cutoff) {
				to_left.push_back(*(sys->get_atom(i)));
				num_to_left++;
			}
			if  (pbc(sys->get_atom(i)->pos, box)[PARALLELDIM] > (rank + 1) * box[PARALLELDIM] / nprocs - skin_cutoff) {
				to_right.push_back(*(sys->get_atom(i)));
				num_to_right++;
			}
		}
		
		// KE = sum(i,1/2 *m(i)*v(i)*v(i))
		for (int k = 0; k < NDIM; ++k) {
			kinetic_energy += 0.5*(sys->get_atom(i)->mass*sys->get_atom(i)->vel[k]*sys->get_atom(i)->vel[k]);
		}
	}
	
	if (nprocs > 1) {
		// Send ghost atoms to neighboring processor
		MPI_Isend(&num_to_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req);
		MPI_Irecv(&num_from_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req+1);
		MPI_Isend(&num_to_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+2);
		MPI_Irecv(&num_from_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+3);
		MPI_Waitall (4, req, stat);

		from_left.resize(num_from_left);
		from_right.resize(num_from_right);
		
		MPI_Isend(&to_left[0], num_to_left, MPI_ATOM, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2);
		MPI_Irecv(&from_left[0], num_from_left, MPI_ATOM, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2+1);
		MPI_Isend(&to_right[0], num_to_right, MPI_ATOM, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+2);
		MPI_Irecv(&from_right[0], num_from_right, MPI_ATOM, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+3);
		MPI_Waitall (4, req2, stat2);

		// Calculate interactions between new (ghost) atoms and atoms on processor
		for (int i=0; i!=sys->natoms(); ++i) {
			for (int j=0; j < num_from_left; ++j) {
				try {
					dE = sys->interact[(*(sys->get_atom(i))).sys_index][from_left[j].sys_index].force_energy(sys->get_atom(i), &from_left[j], &box);
				}
				catch (exception& e) {
					flag_error(e.what(), __FILE__, __LINE__);
					return ILLEGAL_VALUE;
				}
				// Update PE only with interactions from atoms on the left so this energy is not doubly counted in the MPI_Allreduce() below
				potential_energy += dE;
			}
			for (int j=0; j < num_from_right; ++j) {
				try {
					dE = sys->interact[(*(sys->get_atom(i))).sys_index][from_right[j].sys_index].force_energy(sys->get_atom(i), &from_right[j], &box);
				}
				catch (exception& e) {
					flag_error(e.what(), __FILE__, __LINE__);
					return ILLEGAL_VALUE;
				}
			}
		}
	}	

	// Keep track of these on all processors (needed for things like thermostats, etc.)
	MPI_Allreduce (&kinetic_energy, &totKE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce (&potential_energy, &totPE, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	sys->set_total_KE(totKE);
	sys->set_total_PE(totPE);
	if (rank == 0) {
		double totE = totPE + totKE; 
		cout << "KE = " << totKE<< ", PE = " << totPE << ", total = " << totE <<  endl;
	}
	
	// Must wait for all forces to finish calculating before continuing
	MPI_Barrier(MPI_COMM_WORLD);
	return SAFE_EXIT;
}
