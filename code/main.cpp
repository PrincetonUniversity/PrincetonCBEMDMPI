#include "integrator.h"
#include "system.h"
#include "CBEMD.h"

inline void force_calc(System *sys, void (*force)(Atom *a1, Atom *a2, const vector <double> *box)) {
	const vector<double> box = sys->box();
	for (int i=0; i!=sys->natoms(); ++i) {
		for (int j=i+1; j!=sys->natoms(); ++j) {
			//(*force)(sys->get_atom(i), sys->get_atom(j), &sys->box());
			(*force)(sys->get_atom(i), sys->get_atom(j), &box);
		}
	}
}

int main(int argc, char** argv) {

	// make new system    
	System sys1;

	// set system box
	vector<double> box(3,5.0);
	sys1.set_box(box);
	
	Interaction interaction1;

	// add atoms
	Atom atom1;
	atom1.pos[0] = 0.0;
	atom1.pos[1] = 0.0;
	atom1.pos[2] = 0.0;
	atom1.vel[0] = 10.;
	atom1.vel[1] = 10.;
	atom1.vel[2] = 10.;
	atom1.mass = 1.0;
	atom1.type = 1;

	atom::Atom atom2;
	atom2.pos[0] = 1.0;
	atom2.pos[1] = 1.0;
	atom2.pos[2] = 1.0;
	atom2.vel[0] = 0.5;
	atom2.vel[1] = 0.5;
	atom2.vel[2] = 0.5;
	atom2.mass = 1.0;
	atom2.type = 1;
	
	Atom atom_array [] = {atom1, atom2};
	sys1.add_atoms(2, atom_array);
	
	double dt = 0.01;
	Verlet integrator1 (dt);
	int nsteps = 100;
	for (int i = 0; i < nsteps; ++i) {
		force_calc(&sys1, force_serial);
		integrator1.step(&sys1);
		cout << sys1.get_atom(0)->pos[0] << " " << sys1.get_atom(0)->pos[1] << " " <<sys1.get_atom(0)->pos[2] << endl;
	}
	cout << "YIPPEE" << endl;
	return 0;

}
