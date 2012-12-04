#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

//! Given the box size and factors of nprocs, checks which combination generates the most cubic domains
int gen_sets (const vector<int>& factors, const double box[], const int level, double& final_diff, vector<int>& final_breakup, int add);

//! Given a number returns the prime factors (does not count 1 as a prime factor)
vector<int> factorize (const int nprocs);

//! Decomposes the box into domains for each processor to handle
int init_domain_decomp (const vector<double> box, const int nprocs, double widths[], vector<int>& final_breakup);

//! Given the co-ordinates of a point, determines within which domain the point lies
int get_processor (const vector<double> pos, const double widths[], const vector<int>& final_breakup);

//! Generates a mapping of processor id to the x, y, z ids of the domain
int gen_domain_info (const double widths[], const vector<int>& final_breakup, int proc_map[][3]);
