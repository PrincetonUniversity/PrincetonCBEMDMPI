#include "../gtest-1.6.0/include/gtest/gtest.h"
#include "domain_decomp.h"
#include "system.h"
#include "integrator.h"
#include "mpiatom.h"
#include "initialize.h"
#include "read_xml.h"

using namespace std;
using namespace sim_system;
using namespace atom;
using namespace integrator;


class TwoBodyTest : public ::testing::Test {
protected:
    /* Generates a 2 atom system */
    virtual void SetUp () {
	vector<double> box (3, 10.0);
	sys.set_box(box);
        Atom atom0;
        atom0.pos[0] = 2.0;
        atom0.pos[1] = 3.0;
        atom0.pos[2] = 4.0;
        atom0.vel[0] = 1.0;
        atom0.vel[1] = 2.0;
        atom0.vel[2] = 3.0;
	atom0.force[0] = 0.0;
	atom0.force[1] = 0.0;
	atom0.force[2] = 0.0;
        atom0.mass = 1.0;
        atom0.type = 1;
	atom::Atom atom1;
        atom1.pos[0] = 1.0;
        atom1.pos[1] = 2.0;
        atom1.pos[2] = 3.0;
        atom1.vel[0] = 0.1;
        atom1.vel[1] = 0.2;
        atom1.vel[2] = 0.3;
	atom1.force[0] = 0.0;
	atom1.force[1] = 0.0;
	atom1.force[2] = 0.0;
        atom1.mass = 1.0;
        atom1.type = 1;
        Atom atom_array [] = {atom0, atom1};
        sys.add_atoms(2, atom_array);
	sys.set_num_atoms(2);
	sys.final_proc_breakup.reserve(3);
	for (int i=0; i<3; i++) {
	    sys.final_proc_breakup[i] = i+3;
	}
	sys.proc_widths[0] = 10.0/3;
	sys.proc_widths[1] = 2.5;
	sys.proc_widths[2] = 2.0;
    }

    System sys;
};

class ManyBodyTest : public ::testing::Test {
protected:
    /* Generates a many atom system */
    virtual void SetUp () {
	double box_dims[3]={3.0, 4.0, 5.0};
	vector<double> box(3);
	box.assign(box_dims, box_dims+3);
	sys.set_box(box);
	// Add a bunch of random atoms
	for (int i=0; i<3; i++) {
	    for (int j=0; j<4; j++) {
		for (int k=0; k<5; k++) {
		    Atom* atom_ptr = new Atom;
		    int vals[3]={i, j, k};
		    for (int m=0; m<3; m++) {
			atom_ptr->pos[m] = (vals[m]+0.5)*box[m]/6.0;
			atom_ptr->vel[m] = 0.0;
			atom_ptr->force[m] = 0.0;
			atom_ptr->mass = 1.0;
			atom_ptr->type = 1.0;
		    }
		    int* temp = sys.add_atoms(1, atom_ptr);
		    delete atom_ptr;
		}
	    }
	}
	sys.final_proc_breakup.reserve(3);
	for (int i=0; i<3; i++) {
	    sys.proc_widths[i] = 1.0;
	    sys.final_proc_breakup[i] = i+3.0;
	}
	sys.set_num_atoms(3*4*5);
	sys.set_rank(7);
    }

    System sys;
};

TEST_F (ManyBodyTest, SendTableTest) {
    int status=3;
    status = sys.gen_domain_info ();
    ASSERT_EQ (0.0, status);
    EXPECT_EQ (1.0, sys.xyz_id[0]);
    EXPECT_EQ (2.0, sys.xyz_id[1]);
    EXPECT_EQ (0.0, sys.xyz_id[2]);
    status = gen_send_table (&sys);
    ASSERT_EQ (0.0, status);
}

TEST_F (TwoBodyTest, VerletNoForces) {
    Verlet integ_obj (0.1);
    int status;
    for (int i=0; i<10; i++) {
	status = integ_obj.step(&sys);
    }
    ASSERT_EQ (0.0, status);
    EXPECT_DOUBLE_EQ (3.0, sys.get_atom(0)->pos[0]);
    EXPECT_DOUBLE_EQ (5.0, sys.get_atom(0)->pos[1]);
    EXPECT_DOUBLE_EQ (7.0, sys.get_atom(0)->pos[2]);
    EXPECT_DOUBLE_EQ (1.1, sys.get_atom(1)->pos[0]);
    EXPECT_DOUBLE_EQ (2.2, sys.get_atom(1)->pos[1]);
    EXPECT_DOUBLE_EQ (3.3, sys.get_atom(1)->pos[2]);
}

TEST_F (TwoBodyTest, GetProcessorTest) {
    vector<double> pos;
    int proc;
    for (int i=1; i<=3; i++) {
	for (int j=1; j<=4; j++) {
	    for (int k=1; k<=5; k++) {
		pos.clear();
		pos.push_back(1+(i-1)*3.33);
		pos.push_back(1+(j-1)*2.5);
		pos.push_back(1+(k-1)*2);
		proc = get_processor (pos, sys.proc_widths, sys.final_proc_breakup);
		EXPECT_EQ ((i-1)+3*(j-1)+12*(k-1), proc);
	    }
	}
    }
}

TEST_F (TwoBodyTest, VerletNoInitVelocity) {
    Verlet integ_obj (0.1);
    for (int i=0; i<3; i++) {
	sys.get_atom(0)->force[i] = i;
	sys.get_atom(1)->force[i] = i*10;
	sys.get_atom(0)->vel[i] = 0.0;
	sys.get_atom(1)->vel[i] = 0.0;
    }
    sys.get_atom(1)->mass = 2.0;
    integ_obj.step(&sys);
    EXPECT_DOUBLE_EQ (2.0, sys.get_atom(0)->pos[0]);
    EXPECT_DOUBLE_EQ (3.005, sys.get_atom(0)->pos[1]);
    EXPECT_DOUBLE_EQ (4.01, sys.get_atom(0)->pos[2]);
    EXPECT_DOUBLE_EQ (1.0, sys.get_atom(1)->pos[0]);
    EXPECT_DOUBLE_EQ (2.025, sys.get_atom(1)->pos[1]);
    EXPECT_DOUBLE_EQ (3.05, sys.get_atom(1)->pos[2]);

    integ_obj.step(&sys);
    EXPECT_DOUBLE_EQ (2.0, sys.get_atom(0)->pos[0]);
    EXPECT_DOUBLE_EQ (3.02, sys.get_atom(0)->pos[1]);
    EXPECT_DOUBLE_EQ (4.04, sys.get_atom(0)->pos[2]);
    EXPECT_DOUBLE_EQ (1.0, sys.get_atom(1)->pos[0]);
    EXPECT_DOUBLE_EQ (2.1, sys.get_atom(1)->pos[1]);
    EXPECT_DOUBLE_EQ (3.2, sys.get_atom(1)->pos[2]);
}

TEST (DomainDecompTest, OneProcessorTest) {
    vector<double> box;
    box.push_back(10.0);
    box.push_back(10.0);
    box.push_back(10.0);
    double widths[3];
    int nprocs=1, status;
    vector<int> final_breakup;
    status = init_domain_decomp (box, nprocs, widths, final_breakup);
    ASSERT_EQ (0, status);    
    EXPECT_EQ (10.0, widths[0]);
    EXPECT_EQ (10.0, widths[1]);
    EXPECT_EQ (10.0, widths[2]);
    EXPECT_EQ (1, final_breakup[0]);
    EXPECT_EQ (1, final_breakup[1]);
    EXPECT_EQ (1, final_breakup[2]);
}

TEST (DomainDecompTest, RptFactor) {
    vector<double> box;
    box.push_back(3.0);
    box.push_back(4.0);
    box.push_back(2.0);
    double widths[3];
    int nprocs=24, status;
    vector<int> final_breakup;
    status = init_domain_decomp (box, nprocs, widths, final_breakup);
    ASSERT_EQ (0, status);    
    EXPECT_EQ (1.0, widths[0]);
    EXPECT_EQ (1.0, widths[1]);
    EXPECT_EQ (1.0, widths[2]);
    EXPECT_EQ (3, final_breakup[0]);
    EXPECT_EQ (4, final_breakup[1]);
    EXPECT_EQ (2, final_breakup[2]);
}

TEST (DomainDecompTest, UniqueFactor) {
    vector<double> box;
    box.push_back(3.0);
    box.push_back(5.0);
    box.push_back(2.0);
    double widths[3];
    int nprocs=30, status;
    vector<int> final_breakup;
    status = init_domain_decomp (box, nprocs, widths, final_breakup);
    ASSERT_EQ (0, status);    
    EXPECT_EQ (1.0, widths[0]);
    EXPECT_EQ (1.0, widths[1]);
    EXPECT_EQ (1.0, widths[2]);
    EXPECT_EQ (3, final_breakup[0]);
    EXPECT_EQ (5, final_breakup[1]);
    EXPECT_EQ (2, final_breakup[2]);
}

TEST (DomainDecompTest, TwoFactor) {
    vector<double> box;
    box.push_back(3.0);
    box.push_back(1.0);
    box.push_back(2.0);
    double widths[3];
    int nprocs=6, status;
    vector<int> final_breakup;
    status = init_domain_decomp (box, nprocs, widths, final_breakup);
    ASSERT_EQ (0, status);    
    EXPECT_EQ (1.0, widths[0]);
    EXPECT_EQ (1.0, widths[1]);
    EXPECT_EQ (1.0, widths[2]);
    EXPECT_EQ (3, final_breakup[0]);
    EXPECT_EQ (1, final_breakup[1]);
    EXPECT_EQ (2, final_breakup[2]);
}

TEST (DomainDecompTest, OneFactor) {
    vector<double> box;
    box.push_back(1.0);
    box.push_back(5.0);
    box.push_back(1.0);
    double widths[3];
    int nprocs=5, status;
    vector<int> final_breakup;
    status = init_domain_decomp (box, nprocs, widths, final_breakup);
    ASSERT_EQ (0, status);    
    EXPECT_EQ (1.0, widths[0]);
    EXPECT_EQ (1.0, widths[1]);
    EXPECT_EQ (1.0, widths[2]);
    EXPECT_EQ (1, final_breakup[0]);
    EXPECT_EQ (5, final_breakup[1]);
    EXPECT_EQ (1, final_breakup[2]);
}

TEST (DomainDecompTest, TwoFactorLongBox) {
    vector<double> box;
    box.push_back(1.0);
    box.push_back(14.0);
    box.push_back(1.0);
    double widths[3];
    int nprocs=14, status;
    vector<int> final_breakup;
    status = init_domain_decomp (box, nprocs, widths, final_breakup);
    ASSERT_EQ (0, status);    
    EXPECT_EQ (1.0, widths[0]);
    EXPECT_EQ (1.0, widths[1]);
    EXPECT_EQ (1.0, widths[2]);
    EXPECT_EQ (1, final_breakup[0]);
    EXPECT_EQ (14, final_breakup[1]);
    EXPECT_EQ (1, final_breakup[2]);
}

TEST (DomainDecompTest, ThreeFactorLongBox1) {
    vector<double> box;
    box.push_back(1.0);
    box.push_back(10.0);
    box.push_back(3.0);
    double widths[3];
    int nprocs=30, status;
    vector<int> final_breakup;
    status = init_domain_decomp (box, nprocs, widths, final_breakup);
    ASSERT_EQ (0, status);    
    EXPECT_EQ (1.0, widths[0]);
    EXPECT_EQ (1.0, widths[1]);
    EXPECT_EQ (1.0, widths[2]);
    EXPECT_EQ (1, final_breakup[0]);
    EXPECT_EQ (10, final_breakup[1]);
    EXPECT_EQ (3, final_breakup[2]);
}

TEST (DomainDecompTest, ThreeFactorLongBox2) {
    vector<double> box;
    box.push_back(1.0);
    box.push_back(30.0);
    box.push_back(1.0);
    double widths[3];
    int nprocs=30, status;
    vector<int> final_breakup;
    status = init_domain_decomp (box, nprocs, widths, final_breakup);
    ASSERT_EQ (0, status);    
    EXPECT_EQ (1.0, widths[0]);
    EXPECT_EQ (1.0, widths[1]);
    EXPECT_EQ (1.0, widths[2]);
    EXPECT_EQ (1, final_breakup[0]);
    EXPECT_EQ (30, final_breakup[1]);
    EXPECT_EQ (1, final_breakup[2]);
}

TEST (DomainDecompTest, ThreeFactorLongBox3) {
    vector<double> box;
    box.push_back(10.0);
    box.push_back(1.0);
    box.push_back(3.0);
    double widths[3];
    int nprocs=30, status;
    vector<int> final_breakup;
    status = init_domain_decomp (box, nprocs, widths, final_breakup);
    ASSERT_EQ (0, status);    
    EXPECT_EQ (1.0, widths[0]);
    EXPECT_EQ (1.0, widths[1]);
    EXPECT_EQ (1.0, widths[2]);
    EXPECT_EQ (10, final_breakup[0]);
    EXPECT_EQ (1, final_breakup[1]);
    EXPECT_EQ (3, final_breakup[2]);    
}
			   			   


// the following use mpi in the tests
TEST (ReadXMLTest, AtomPositions) {
    int argc = 5;
    char *argv[] = {"dummy"};
    int check = start_mpi (argc, argv);
    System sys1;
    int status;
    status = initialize_from_files ("sample.xml", "energy_params.dat",&sys1);
    //status = read_xml("sample.xml", 4, 0, &sys1);
    //cout << status << endl;
    cout << "Total Atoms " << sys1.total_atoms() << " " << sys1.natoms() << endl;
    cout << "x,y,z coords of atom 1 " << sys1.get_atom(0)->pos[0] << " " << sys1.get_atom(0)->pos[1] << " " << sys1.get_atom(0)->pos[2] << endl;
    ASSERT_EQ(5,sys1.get_atom(0)->pos[0]);
}

TEST (ReadXMLTest, AtomVelocities) { 
    int argc = 5;
    char *argv[] = {"dummy"};
    //int check = start_mpi (argc, argv);
    System sys1;
    int status;
    status = initialize_from_files ("sample.xml", "energy_params.dat",&sys1);
    //status = read_xml("sample.xml", 4, 0, &sys1);
    //cout << status << endl;
    cout << "Total Atoms " << sys1.total_atoms() << " " << sys1.natoms() << endl;
    cout << "x,y,z, velocities of atom 1 " << sys1.get_atom(0)->vel[0] << " " << sys1.get_atom(0)->vel[1] << " " << sys1.get_atom(0)->vel[2] << endl;
    ASSERT_EQ(0,sys1.get_atom(0)->vel[0]);
}

TEST (ReadXMLTest, BoxVolume) {
    int argc = 5;
    char *argv[] = {"dummy"};
    //int check = start_mpi (argc, argv);
    System sys1;
    int status;
    status = initialize_from_files ("sample.xml", "energy_params.dat",&sys1);
    //status = read_xml("sample.xml", 4, 0, &sys1);
    //cout << status << endl;
    for (int i = 0; i < 3; i++) {
        cout<<sys1.box()[i] << " " ;
    }
    cout << endl;
    ASSERT_EQ(1000,sys1.box()[0]*sys1.box()[1]*sys1.box()[2]);
    MPI_Finalize();
}

// probably check that each error works

int main (int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
