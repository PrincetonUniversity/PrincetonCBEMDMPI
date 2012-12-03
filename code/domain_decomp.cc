#include "domain_decomp.h"

//!< Given the box size and factors of nprocs, checks which combination generates the most cubic domains
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
} //gen_sets ends

//!< Given a number returns the prime factors (does not count 1 as a prime factor)
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

//!< Decomposes the box into domains for each processor to handle
int init_domain_decomp (const vector<double> box, const int nprocs, double widths[], vector<int>& final_breakup) {

    int status;
    double box_dims[3];
    double final_diff=10000;
    vector<int> factors;

    final_breakup.reserve(3);
    for (vector<int>::size_type it=0; it!=final_breakup.capacity(); it++) {
	final_breakup.push_back(-1);
    }
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
} //init_domain_decomp ends

//!< Given the co-ordinates of a point, determines within which domain the point lies
int get_processor (const vector<double> pos, const double widths[], const vector<int>& final_breakup) {

    int domain_id, x_id, y_id, z_id;
    x_id = floor(pos.at(0)/widths[0]);
    y_id = floor(pos.at(1)/widths[1]);
    z_id = floor(pos.at(2)/widths[2]);
    domain_id = x_id + y_id*final_breakup[0] + z_id*final_breakup[0]*final_breakup[1];
    return domain_id;
} //get_processor ends

//!< Generates a mapping of processor id to the x, y, z ids of the domain
int gen_domain_info (const double widths[], const vector<int>& final_breakup, int proc_map[][3]) {

    int domain_id;
    for (int x_id=0; x_id<final_breakup[0]; x_id++) {
	for (int y_id=0; y_id<final_breakup[1]; y_id++) {
	    for (int z_id=0; z_id<final_breakup[2]; z_id++) {
		domain_id = x_id + y_id*final_breakup[0] + z_id*final_breakup[0]*final_breakup[1];
		proc_map[domain_id][0] = x_id;
		proc_map[domain_id][1] = y_id;
		proc_map[domain_id][2] = z_id;
	    }
	}
    }
    return 0;
} // gen_domain_info ends
