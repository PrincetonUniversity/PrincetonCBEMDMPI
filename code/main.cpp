/*!
 \brief Test MD driver program for serial implementation of code
 \file main.cpp
 \authors{Nathan A. Mahynski, Carmeline Dsilva, Arun L. Prabhu, George Khoury}
 **/

#include "integrator.h"
#include "system.h"
#include "CBEMD.h"

/*!
 \brief Function to perform the force calculation
 \param [in] a1 Pointer to Atom 1 involved in force calculation
 \param [in] a2 Pointer to Atom 2 involved in force calculation
 \param [in] box Pointer to the system box
*/
inline void force_calc(System *sys, void (*force)(Atom *a1, Atom *a2, const vector <double> *box)) {
	const vector<double> box = sys->box();
	for (int i=0; i!=sys->natoms(); ++i) {
		for (int j=i+1; j!=sys->natoms(); ++j) {
			//(*force)(sys->get_atom(i), sys->get_atom(j), &sys->box());
			(*force)(sys->get_atom(i), sys->get_atom(j), &box);
			//cout << sys->natoms() << endl;
			//cout << sys->get_atom(0)->force[0] << " " << sys->get_atom(0)->force[1] << " " <<sys->get_atom(0)->force[2] << endl;
			
		}
	}
}

/*!
 \brief Main setting up a few atoms and their velocities and performing MD 
*/
int main(int argc, char** argv) {

	// make new system    
	System sys1;

	// set system box
	vector<double> box(3,10.0);
	sys1.set_box(box);
	
	Interaction interaction1;

	// add atoms
	Atom atom1;
	atom1.pos[0] = 2.0;
	atom1.pos[1] = 3.0;
	atom1.pos[2] = 4.0;
	atom1.vel[0] = 1.0;
	atom1.vel[1] = 2.0;
	atom1.vel[2] = 3.0;
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
	int nsteps = 10;
	for (int i = 0; i < nsteps; ++i) {
		force_calc(&sys1, force_serial);
		integrator1.step(&sys1);
		cout << sys1.get_atom(0)->pos[0] << " " << sys1.get_atom(0)->pos[1] << " " <<sys1.get_atom(0)->pos[2] << endl;
	}
	cout << "YIPPEE" << endl;  // Happiness if successful.
	return 0;

}
