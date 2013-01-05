#include "../gtest-1.6.0/include/gtest/gtest.h"
#include "domain_decomp.h"
#include "system.h"
#include "integrator.h"
#include "mpiatom.h"
#include "gtest_mpi.h"

using namespace std;
using namespace sim_system;
using namespace atom;
using namespace integrator;


int main (int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    testing::TestEventListeners& listeners = testing::UnitTest::GetInstance()->listeners();
    listeners.Append(new MPI_PrettyUnitTestResultPrinter);
    return RUN_ALL_TESTS();
}
