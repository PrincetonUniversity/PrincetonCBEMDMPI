#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>
#include "system.h"

using namespace std;
using namespace sim_system;

//! Given the box size and factors of nprocs, checks which combination generates the most cubic domains
int gen_sets (const vector<int>& factors, const double box[], const int level, double& final_diff, vector<int>& final_breakup, int add);

//! Given a number returns the prime factors (does not count 1 as a prime factor)
vector<int> factorize (const int nprocs);

//! Decomposes the box into domains for each processor to handle
int init_domain_decomp (const vector<double> box, const int nprocs, double widths[], vector<int>& final_breakup);

//! Given the co-ordinates of a point, determines within which domain the point lies
inline int get_processor (const vector<double> pos, const double widths[], const vector<int>& final_breakup);

//! Given the x, y, z ids of a domain, determines the domain id (useful for locating neighbouring domains)
inline int get_processor (const int x_id, const int y_id, const int z_id, const vector<int>& final_breakup);

//! Generates the x,y,z ids for each processor and the absolute extents of the domain
int gen_domain_info (System *sys, const double widths[], const vector<int>& final_breakup, const int rank);

//! Given a domain_id specifies the x, y, z ids of the domain
int get_xyz_ids (const int domain_id, const vector<int>& final_breakup, int xyz_id[]);

//! Generates the lists of molecules that need to be passed to other processors
int gen_send_lists (System *sys, const int rank, const double skin_cutoff);

//! Given the rank, generates the list of its neighbours
int gen_send_table (System *sys);

//! Generates the skin cutoff distance based on the largest cutoff in the system
double gen_skin_cutoff ( ... /*Need a list of all the cutoffs for the different interactions*/ );

//! Generates the list of neighbours a particle should be sent to based on the borders its near
void gen_goes_to (const int is_near_border[], vector<int>& goes_to, const int ndims);

//! Computes the exponentiation of an integer by an integral power
int power (int base, int exponent);
