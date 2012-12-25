/**
 \file test_all.cpp
 \authors {Arun L. Prabhu, Carmeline Dsilva, Nathan Mahynski, George Khoury}
 **/
 
#include <gtest/gtest.h>
#include "atom.h"
#include "integrator.h"
#include "system.h"

TEST(IntegratorTest, Test1) {
       integrator::Verlet test1 (0.1);
}

TEST(AtomTest, Test1) {
       atom::Atom atom1;
       atom1.pos[0] = 0.0;
       atom1.pos[1] = 0.0;
       atom1.pos[2] = 0.0;
       
       atom1.vel[0] = 0.5;
       atom1.vel[1] = 0.5;
       atom1.vel[2] = 0.5;
       
       atom1.mass = 1.0;
       atom1.type = 1;
       
       EXPECT_EQ(atom1.pos[0], 0.0);

}

TEST(SystemTest, Test1) {
       sim_system::System sys1;
       atom::Atom atom1;
       sys1.add_atoms(1, &atom1);
       
       atom::Atom atom2;
       sys1.add_atoms(1, &atom2);
       
       EXPECT_EQ(sys1.natoms(), 2);
       
       std::vector<int> to_delete (1, 0);
       sys1.delete_atoms(to_delete);
       
       EXPECT_EQ(sys1.natoms(), 1);
       
}

int main(int argc, char** argv) {
       testing::InitGoogleTest(&argc, argv);
       return RUN_ALL_TESTS();
}
