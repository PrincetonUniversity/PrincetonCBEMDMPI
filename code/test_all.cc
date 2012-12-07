#include "../gtest-1.6.0/include/gtest/gtest.h"
#include "domain_decomp.h"
#include "system.h"
#include "integrator.h"

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
    }

    System sys;
};

TEST_F (TwoBodyTest, ForceSerial) {
    vector<double> box=sys.box();
    /* Moves particle 0 within the cutoff distance of particle 1 so that a force will be calculate */
    sys.get_atom(0)->pos[0] = 3.0;
    sys.get_atom(0)->pos[1] = 2.0;
    sys.get_atom(0)->pos[2] = 3.0;
    force_serial(sys.get_atom(0), sys.get_atom(1), &box);
    /* Checks if the forces are correct if the particles are within the cutoff */
    EXPECT_FLOAT_EQ (-2.3999952e-6, sys.get_atom(0)->force[0]);
    EXPECT_FLOAT_EQ (0.0, sys.get_atom(0)->force[1]);
    EXPECT_FLOAT_EQ (0.0, sys.get_atom(0)->force[2]);
    EXPECT_FLOAT_EQ (2.3999952e-6, sys.get_atom(1)->force[0]);
    EXPECT_FLOAT_EQ (0.0, sys.get_atom(1)->force[1]);
    EXPECT_FLOAT_EQ (0.0, sys.get_atom(1)->force[2]);
    /* Resets the forces to zero for the next part of the test */
    sys.get_atom(0)->force[0] = 0.0;
    sys.get_atom(1)->force[0] = 0.0;
    /* Moves particle 0 outside the cutoff distance from particle 1 so that a force should not be calculated */
    sys.get_atom(0)->pos[0] = 5.0;
    sys.get_atom(0)->pos[1] = 2.0;
    sys.get_atom(0)->pos[2] = 3.0;
    force_serial(sys.get_atom(0), sys.get_atom(1), &box);
    /* Checks if the forces are correct if the particles are outside the cutoff */
    EXPECT_DOUBLE_EQ (0.0, sys.get_atom(0)->force[0]);
    EXPECT_DOUBLE_EQ (0.0, sys.get_atom(0)->force[1]);
    EXPECT_DOUBLE_EQ (0.0, sys.get_atom(0)->force[2]);
    EXPECT_DOUBLE_EQ (0.0, sys.get_atom(1)->force[0]);
    EXPECT_DOUBLE_EQ (0.0, sys.get_atom(1)->force[1]);
    EXPECT_DOUBLE_EQ (0.0, sys.get_atom(1)->force[2]);
}

TEST_F (TwoBodyTest, VerletNoForces) {
    Verlet integ_obj (0.1);
    for (int i=0; i<10; i++) {
	integ_obj.step(&sys);
    }
    EXPECT_DOUBLE_EQ (3.0, sys.get_atom(0)->pos[0]);
    EXPECT_DOUBLE_EQ (5.0, sys.get_atom(0)->pos[1]);
    EXPECT_DOUBLE_EQ (7.0, sys.get_atom(0)->pos[2]);
    EXPECT_DOUBLE_EQ (1.1, sys.get_atom(1)->pos[0]);
    EXPECT_DOUBLE_EQ (2.2, sys.get_atom(1)->pos[1]);
    EXPECT_DOUBLE_EQ (3.3, sys.get_atom(1)->pos[2]);
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

int main (int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
