#include "../gtest-1.6.0/include/gtest/gtest.h"
#include "domain_decomp.h"

TEST (DomainDecompTest, Test1) {
    vector<double> box;
    box.push_back(3.0);
    box.push_back(4.0);
    box.push_back(2.0);
    double widths[3];
    int nprocs=24, status;
    vector<int> final_breakup;
    status = init_domain_decomp (box, nprocs, widths, final_breakup);
    ASSERT_EQ (1.0, widths[0]);
    ASSERT_EQ (1.0, widths[1]);
    ASSERT_EQ (1.0, widths[2]);
}

int main (int argc, char* argv[]) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
