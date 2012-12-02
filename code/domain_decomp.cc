#include "domain_decomp.h"

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

vector <int> factorize (const int nprocs) {
    vector <int> factors;
    int unfactored=nprocs;
    factors.push_back(1);
    factors.push_back(1);
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

int init_domain_decomp (const vector<double> box, const int nprocs, double widths[]) {

    int status;
    double box_dims[3];
    double final_diff=10000;
    vector<int> final_breakup;
    vector<int> factors;

    final_breakup.reserve(3);
    for (vector<int>::size_type it=0; it!=final_breakup.capacity(); it++) {
	final_breakup.push_back(-1);
    }
    for (int i=0; i<3; i++) {
	box_dims[i] = box[i];
    }

    factors = factorize (nprocs);
    status = gen_sets(factors, box_dims, 0, final_diff, final_breakup, 0);

    for (int i=0; i<3; i++) {
	widths[i] = box_dims[i] / final_breakup[i];
    }

    return 0;
} //init_domain_decomp ends

int get_processor () {

    return 0;
} //get_processor ends
