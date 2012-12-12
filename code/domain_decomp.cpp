#include "domain_decomp.h"
#include "mpiatom.h"

//! Given the box size and factors of nprocs, checks which combination generates the most cubic domains
int gen_sets (const vector<int>& factors, const double box[], const int level, double& final_diff, vector<int>& final_breakup, int add) {

    int status;
    vector<int> remaining;
    double widths[3];
    static int breakup[3]={1,1,1};

    if (factors.size() == 0) {
	for (int m=0; m<3; m++) {
	    widths[m] = box[m] / breakup[m];
	}
	double diff = (widths[1]-widths[0])*(widths[1]-widths[0]) + (widths[2]-widths[1])*(widths[2]-widths[1]) + (widths[0]-widths[2])*(widths[0]-widths[2]);
	if (diff <= final_diff) {
	    for (vector<int>::size_type m=0; m!=final_breakup.capacity(); m++) {
		final_breakup.at(m) = breakup[m];
	    }
	    final_diff = diff;
	}
	return 0;
    } else if (level != 0 && add != 0) {
	for (int j=0; j<factors.size(); j++) {
	    breakup[level-1] *= factors[j];
	    remaining = factors;
	    remaining.erase (remaining.begin()+j);
	    status = gen_sets(remaining, box, level, final_diff, final_breakup, add-1);
	    breakup[level-1] /= factors[j];
	}
    } else if (level == 0 && add == 0) {
	for (int m=0; m<3; m++) {
	    breakup[m] = 1;
	}
	for (int i=1; i<=(factors.size()-2); i++) {
	    status = gen_sets(factors, box, 1, final_diff, final_breakup, i);
	}
    } else if (level == 2 && add == 0) {
	status = gen_sets(factors, box, level+1, final_diff, final_breakup, factors.size());
    } else if (level != 0 && add == 0) {
	for (int k=1; k<=(factors.size()+level-2); k++) {
	    status = gen_sets(factors, box, level+1, final_diff, final_breakup, k);
	}
    }
    return 0;
} //gen_sets ends

//! Given a number returns the prime factors (does not count 1 as a prime factor)
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
    return factors;
} // factorize ends

//! Decomposes the box into domains for each processor to handle
int init_domain_decomp (const vector<double> box, const int nprocs, double widths[], vector<int>& final_breakup) {

    int status;
    double box_dims[3];
    double final_diff=10000;
    vector<int> factors;

    final_breakup.resize(3, -1);
	
    for (int i=0; i<3; i++) {
	box_dims[i] = box[i];
    }

    factors = factorize (nprocs);
    // Two 1's added to the factors to ensure deompositions like 1,1,6 are found
    factors.push_back(1);
    factors.push_back(1);

    // gen_sets is a recursive function, initial call must have level=0 and add=0
    status = gen_sets(factors, box_dims, 0, final_diff, final_breakup, 0);
    for (int i=0; i<3; i++) {
	widths[i] = box_dims[i] / final_breakup[i];
    }

    return 0;
} // init_domain_decomp ends

//! Given the co-ordinates of a point, determines within which domain the point lies
int get_processor (const vector<double> pos, const double widths[], const vector<int>& final_breakup) {

    int domain_id, x_id, y_id, z_id;
    x_id = floor(pos.at(0)/widths[0]);
    y_id = floor(pos.at(1)/widths[1]);
    z_id = floor(pos.at(2)/widths[2]);
    domain_id = x_id + y_id*final_breakup[0] + z_id*final_breakup[0]*final_breakup[1];
    return domain_id;
}

int get_processor (const int x_id, const int y_id, const int z_id, const vector<int>& final_breakup) {

    int domain_id = x_id + y_id*final_breakup[0] + z_id*final_breakup[0]*final_breakup[1];
    return domain_id;
}

//! Generates the x,y,z ids for each processor and the absolute extents of the domain
int gen_domain_info (System *sys, const double widths[], const vector<int>& final_breakup, const int rank) {

    int domain_id;
    for (int x_id=0; x_id<final_breakup[0]; x_id++) {
	for (int y_id=0; y_id<final_breakup[1]; y_id++) {
	    for (int z_id=0; z_id<final_breakup[2]; z_id++) {
		domain_id = x_id + y_id*final_breakup[0] + z_id*final_breakup[0]*final_breakup[1];
		if (domain_id == rank) {
		    sys->xyz_id[0] = x_id;
		    sys->xyz_id[1] = y_id;
		    sys->xyz_id[2] = z_id;
		    sys->xyz_limits[0][0] = x_id*widths[0];
		    sys->xyz_limits[0][1] = (x_id+1)*widths[0];
		    sys->xyz_limits[1][0] = y_id*widths[1];
		    sys->xyz_limits[1][1] = (y_id+1)*widths[1];
		    sys->xyz_limits[2][0] = z_id*widths[2];
		    sys->xyz_limits[2][1] = (z_id+1)*widths[2];
		    return 0;
		}
	    }
	}
    }
    return 1;
} // gen_domain_info ends

//! Given a domain_id specifies the x, y, z ids of the domain
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

//! Generates the lists of molecules that need to be passed to other processors
int gen_send_lists (System *sys, const int rank, const double skin_cutoff) {

    const int ndims=3;
    /* Since a particle can have only 3 realtionships to a dimension of the box, we define
       0 = in the middle (itm),
       1 = near lower bound (nlb)
       2 = near upper bound (nub) */
    const int itm=0, nlb=1, nub=2;
    int is_near_border[ndims];
    vector<int> goes_to;

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
	    /* Add current atom to those rows given by the values of goes_to
	       send_lists would have to be a vecor of vectors in class System with first dimension 26
	    */
	}
    }
    return 0;
}

//! Given the rank, generates the list of its neighbours
/* Assumes 3D system, for different number of dimensions need to make this a recursive function */
int gen_send_table (System *sys) {

    const int nvals=3;
    int xyz_id[3], domain_id, value;
    for (int i=0; i<3; i++) {
	xyz_id[i] = sys->xyz_id[i];
    }
    for (int i=-1; i<=1; i++) {
	for (int j=-1; j<=1; j++) {
	    for (int k=-1; k<=1; k++) {
		/* Need to include final_breakup and check that
		   if we go greater than final_breakup or less than 0 to wrap around */
		domain_id = i+xyz_id[0] + (j+xyz_id[1])*sys->final_proc_breakup[0] + (k+xyz_id[2])*sys->final_proc_breakup[0]*sys->final_proc_breakup[1];
		if (i*j*k != 0) {
		    value = (1.5*abs(i)+0.5*i)+(1.5*abs(j)+0.5*j)*nvals+(1.5*abs(k)+0.5*k)*nvals*nvals - 1;
		    sys->send_table[value] = domain_id;
		}
	    }
	}
    }
    return 0;
}

//! Generates the skin cutoff distance based on the largest cutoff in the system
double gen_skin_cutoff (/*Need a list of all the cutoffs for the different interactions*/ ) {
    double skin_cutoff;

    return skin_cutoff;
}

//! Generates the list of neighbours a particle should be sent to based on the borders its near
void gen_goes_to (const int is_near_border[], vector<int>& goes_to, const int ndims) {

    int remaining[ndims], value;
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

//! Computes the exponentiation of an integer by an integral power
int power (int base, int exponent) {
    int result = 1;
    for (int i=0; i<exponent; i++) {
	result *= base;
    }
    return result;
}
