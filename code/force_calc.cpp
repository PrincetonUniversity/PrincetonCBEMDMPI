#include "force_calc.h"

using namespace std;
using namespace sim_system;

int force_calc(System *sys) { // pass interaction array
	const vector<double> box = sys->box();
	for (int i=0; i!=sys->total_atoms(); ++i) {
		for (int j=i+1; j!=sys->total_atoms(); ++j) {
		    if ( (*(sys->get_atom(i))).sys_index != (*(sys->get_atom(j))).sys_index ) {
			sys->interact[(*(sys->get_atom(i))).sys_index][(*(sys->get_atom(j))).sys_index].force_energy(sys->get_atom(i), sys->get_atom(j), &box);
		    }
		}
	}
	return 0;
}
