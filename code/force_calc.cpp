#include "system.h"

using namespace std;
using namespace sim_system;

inline void force_calc(System *sys, Interaction *interaction) { // pass interaction array
	const vector<double> box = sys->box();
	for (int i=0; i!=sys->natoms(); ++i) {
		for (int j=i+1; j!=sys->natoms(); ++j) {
			//(*force)(sys->get_atom(i), sys->get_atom(j), &sys->box());
			interaction->force_energy(sys->get_atom(i), sys->get_atom(j), &box);
			//cout << sys->natoms() << endl;
			//cout << sys->get_atom(0)->force[0] << " " << sys->get_atom(0)->force[1] << " " <<sys->get_atom(0)->force[2] << endl;
			
		}
	}
}
