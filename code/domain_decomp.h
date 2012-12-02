#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <vector>

using namespace std;

int gen_sets (const vector<int>& factors, const double box[], const int level, double& final_diff, vector<int>& final_breakup, int add);

vector<int> factorize (const int nprocs);

int init_domain_decomp (const vector<double> box, const int nprocs, double widths[]);

int get_processor ();
