/*!
 \file domain_decomp.cpp
 \brief Source code containing functions for performing domain decomposition
 \author Arun L. Prabhu
*/

#include "domain_decomp.h"

/*!
 Given the box size and factors of nprocs, this recursive function checks which combination generates the most cubic domains
 \param [in] factors WHAT DOES IT DO
 \param [in] box WHAT DOES IT DO
 \param [in] level WHAT DOES IT DO
 \param [in] final_breakup WHAT DOES IT DO
*/
void gen_sets (const vector<int>& factors, const double box[], const int level, double& final_diff, vector<int>& final_breakup, int add) {

    vector<int> remaining;
    double widths[NDIM];
    static int breakup[NDIM]={1,1,1};
    
    if (factors.size() == 0) {
		for (int m=0; m<NDIM; m++) {
			widths[m] = box[m] / breakup[m];
		}
		double diff = (widths[1]-widths[0])*(widths[1]-widths[0]) + (widths[2]-widths[1])*(widths[2]-widths[1]) + (widths[0]-widths[2])*(widths[0]-widths[2]);
		if (diff <= final_diff) {
			for (vector<int>::size_type m=0; m!=final_breakup.capacity(); m++) {
				final_breakup.at(m) = breakup[m];
			}
			final_diff = diff;
		}
		return;
    } else if (level != 0 && add != 0) {
		for (unsigned int j=0; j<factors.size(); j++) {
			breakup[level-1] *= factors[j];
			remaining = factors;
			remaining.erase (remaining.begin()+j);
			gen_sets(remaining, box, level, final_diff, final_breakup, add-1);
			breakup[level-1] /= factors[j];
		}
    } else if (level == 0 && add == 0) {
		for (int m=0; m<NDIM; m++) {
			breakup[m] = 1;
		}
		for (int i=1; i<=((int)factors.size()-2); i++) {
			gen_sets(factors, box, 1, final_diff, final_breakup, i);
		}
    } else if (level == 2 && add == 0) {
		gen_sets(factors, box, level+1, final_diff, final_breakup, factors.size());
    } else if (level != 0 && add == 0) {
		for (int k=1; k<=((int)factors.size()+level-2); k++) {
			gen_sets(factors, box, level+1, final_diff, final_breakup, k);
		}
    }
    return;
} //gen_sets ends

/*!
 Given a number returns the prime factors of that number.
 This function does not count 1 as a prime factor.
 \param [in] nprocs number of processors
*/

vector <int> factorize (const int nprocs) {

    vector <int> factors;
    int unfactored=nprocs;
    int i=2;
    while (i<=unfactored) {
	if ((unfactored%i) == 0) {
	    factors.push_back(i);
	    unfactored /= i;
	    i=2;
	} else {
	    i++;
	}
    }
    if (factors.size() == 0) {
	factors.push_back(1);
    }
    return factors;
} // factorize ends

/*!
 Decomposes the box into domains for each processor to handle
 \param [in] box WHAT DOES IT DO
 \param [in] nprocs number of processors
 \param [in] widths WHAT DOES IT DO
 \param [out] final_breakup WHAT DOES IT DO
*/
int init_domain_decomp (const vector<double> box, const int nprocs, double widths[], vector<int>& final_breakup) {

    double box_dims[NDIM];
    double final_diff=10000;
    vector<int> factors;

    final_breakup.resize(NDIM, -1);
	
    for (int i=0; i<NDIM; i++) {
	box_dims[i] = box[i];
    }

    factors = factorize (nprocs);
    // Two 1's added to the factors to ensure deompositions like 1,1,6 are found
    factors.push_back(1);
    factors.push_back(1);

    // gen_sets is a recursive function, initial call must have level=0 and add=0
    gen_sets(factors, box_dims, 0, final_diff, final_breakup, 0);
    for (int i=0; i<NDIM; i++) {
	widths[i] = box_dims[i] / final_breakup[i];
    }

    return 0;
} // init_domain_decomp ends

/*! Given the coordinates of a point, determines within which domain the point lies. This function is overloaded.
 \param [in] pos the vector containing the coordinates of a point
 \param [in] sys system passed to be able to utilize final domain decomposition
*/
int get_processor (const vector<double> pos, const System *sys) {

	vector <double> inbox = pbc (pos, sys->box());
    int domain_id, x_id, y_id, z_id;
    x_id = floor(inbox[0]/sys->proc_widths[0]);
    y_id = floor(inbox[1]/sys->proc_widths[1]);
    z_id = floor(inbox[2]/sys->proc_widths[2]);
    domain_id = x_id + y_id*sys->final_proc_breakup[0] + z_id*sys->final_proc_breakup[0]*sys->final_proc_breakup[1];
    return domain_id;
}

/*! Given the coordinates of a point, determines within which domain the point lies. This function is overloaded.
 \param [in] x_id WHAT DOES IT DO
 \param [in] y_id WHAT DOES IT DO
 \param [in] z_id WHAT DOES IT DO
 \param [in] fina_breakup WHAT DOES IT DO
*/
int get_processor (const int x_id, const int y_id, const int z_id, const vector<int>& final_breakup) {

    int domain_id = x_id + y_id*final_breakup[0] + z_id*final_breakup[0]*final_breakup[1];
    return domain_id;
}

/*! Given a domain_id specifies the x, y, z ids of the domain
 \param [in] domain_id WHAT DOES IT DO
 \param [in] final_breakup WHAT DOES IT DO
 \param [in] xyz_id WHAT DOES IT DO
*/
int get_xyz_ids (const int domain_id, const vector<int>& final_breakup, int xyz_id[]) {

    int value;
    for (int x_id=0; x_id<final_breakup[0]; x_id++) {
	for (int y_id=0; y_id<final_breakup[1]; y_id++) {
	    for (int z_id=0; z_id<final_breakup[2]; z_id++) {
		value = x_id + y_id*final_breakup[0] + z_id*final_breakup[0]*final_breakup[1];
		if (value == domain_id) {
		    xyz_id[0] = x_id;
		    xyz_id[1] = y_id;
		    xyz_id[2] = z_id;
		    return 0;
		}
	    }
	}
    }
    return 1;
}

/*! Generates the lists of molecules that need to be passed to other processors
 \param [in,out] sys WHAT DOES IT DO
 \param [in] rank rank of processor
 \param [in] skin_cutoff WHAT DOES IT DO
*/
int gen_send_lists (System *sys) {
    const int ndims=NDIM;
    const double skin_cutoff=sys->max_rcut();
    /* Since a particle can have only 3 relationships to a dimension of the box, we define
       0 = in the middle (itm),
       1 = near lower bound (nlb)
       2 = near upper bound (nub) */
    const int itm=0, nlb=1, nub=2;
    vector<int> is_near_border;
    is_near_border.resize(ndims);
    vector<int> goes_to;

    sys->send_lists.clear();
    for (int i=0; i<NNEIGHBORS; i++) {
		sys->send_list_size[i] = 0;
    }
    for (int i=0; i < sys->natoms(); i++) {
		goes_to.clear();
		for (int j=0; j<ndims; j++) {
			if (sys->get_atom(i)->pos[j] < (sys->xyz_limits[j][0]+skin_cutoff)) {
				is_near_border[j] = nlb;
			} else if  (sys->get_atom(i)->pos[j] > (sys->xyz_limits[j][1]-skin_cutoff)) {
				is_near_border[j] = nub;
			} else {
				is_near_border[j] = itm;
			}
		}
		gen_goes_to(is_near_border, goes_to, ndims);
		for (vector<int>::iterator iter=goes_to.begin(); iter!=goes_to.end(); iter++) {
			sys->send_lists[*iter].push_back(*(sys->get_atom(i)));
			sys->send_list_size[*iter]++;
		}
    }
    return 0;
}

/*! Given the rank, generates the list of its neighbours
 Assumes 3D system, for different number of dimensions need to make this a recursive function
 \param [in] sys WHAT DOES IT DO
*/
int gen_send_table (System *sys) {

    const int nvals=NDIM;
    int xyz_id[NDIM], ngh_xyz_id[NDIM], domain_id, value;

    for (int i=0; i<NDIM; i++) {
	xyz_id[i] = sys->xyz_id[i];
    }
    for (int i=-1; i<=1; i++) {
	for (int j=-1; j<=1; j++) {
	    for (int k=-1; k<=1; k++) {
		ngh_xyz_id[0] = i + xyz_id[0];
		ngh_xyz_id[1] = j + xyz_id[1];
		ngh_xyz_id[2] = k + xyz_id[2];
		for (int m=0; m<NDIM; m++) {
		    if (ngh_xyz_id[m] < 0) {
			ngh_xyz_id[m] = sys->final_proc_breakup[m]-1;
		    } else if (ngh_xyz_id[m] == sys->final_proc_breakup[m]) {
			ngh_xyz_id[m] = 0;
		    }
		}
		domain_id = ngh_xyz_id[0] + ngh_xyz_id[1]*sys->final_proc_breakup[0] + ngh_xyz_id[2]*sys->final_proc_breakup[0]*sys->final_proc_breakup[1];
		if (abs(i)+abs(j)+abs(k) != 0) {
		    value = (1.5*abs(i)+0.5*i)+(1.5*abs(j)+0.5*j)*nvals+(1.5*abs(k)+0.5*k)*nvals*nvals - 1;
		    sys->send_table[value] = domain_id;
		}
	    }
	}
    }
    return 0;
}

/*! Generates the skin cutoff distance based on the largest cutoff in the system
*/
double gen_skin_cutoff (/*Need a list of all the cutoffs for the different interactions*/ ) {
    double skin_cutoff=0.5;

    return skin_cutoff;
}

/*! Generates the list of neighbours a particle should be sent to based on the borders its near
 \param [in] is_near_border WHAT DOES IT DO
 \param goes_to [in] WHAT DOES IT DO
 \param [in] ndims WHAT DOES IT DO
*/
void gen_goes_to (const vector<int>& is_near_border, vector<int>& goes_to, const int ndims) {

    int value=0;
    vector<int> remaining;
    remaining.resize(ndims);
    /* Since a particle can have only 3 realtionships to a dimension
       near lower bound, near upper bound or in the middle. Hence the variable nvals */
    const static int nvals=3;
    for (int i=0; i<ndims; i++) {
	remaining[i] = is_near_border[i];
	value += is_near_border[i]*power(nvals, i);
    }
    value--;
    if (value == -1) {
	return;
    } else {
	goes_to.push_back(value);
	for (int j=0; j<ndims; j++) {
	    if (remaining[j] != 0) {
		remaining[j] = 0;
		gen_goes_to (remaining, goes_to, ndims);
		remaining[j] = is_near_border[j];
	    }
	}
	return;
    }
}

/*! Computes the exponentiation of an integer by an integral power
 \param [in] base the base value to be raised to a power
 \param [in] exponent the exponent of the power base is raised to
*/
int power (int base, int exponent) {
    int result = 1;
    for (int i=0; i<exponent; i++) {
	result *= base;
    }
    return result;
}

/*! Communicates the atoms in the skin regions of the processors to the appropriate neighbours
 \param sys [in,out] System to be evaluated
*/
int communicate_skin_atoms (System *sys) {

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Request req[52], req2[52];
    MPI_Status stat[52];
    for (int i=0; i<NNEIGHBORS; i++) {
	MPI_Isend (&sys->send_list_size[i], 1, MPI_INT, sys->send_table[i], sys->rank(), MPI_COMM_WORLD, &req[i]);
	MPI_Irecv (&sys->get_list_size[i], 1, MPI_INT, sys->send_table[i], sys->send_table[i], MPI_COMM_WORLD, &req[NNEIGHBORS+i]);
    }
    MPI_Waitall (52, req, stat);
    
    for(int i=0; i<NNEIGHBORS; i++) {
	sys->get_lists[i].reserve(sys->get_list_size[i]);
    }

    for (int i=0; i<NNEIGHBORS; i++) {
	MPI_Isend (&sys->send_lists[i].front(), sys->send_list_size[i], MPI_ATOM, sys->send_table[i], sys->rank(), MPI_COMM_WORLD, &req2[i]);
	MPI_Irecv (&sys->get_lists[i].front(), sys->get_list_size[i], MPI_ATOM, sys->send_table[i], sys->send_table[i], MPI_COMM_WORLD, &req2[NNEIGHBORS+i]);
    }
    MPI_Waitall (52, req2, stat);

    for (int i=0; i<NNEIGHBORS; i++) {
	sys->add_ghost_atoms(sys->get_list_size[i], &(sys->get_lists[i].front()));
    }
    return 0;
}
