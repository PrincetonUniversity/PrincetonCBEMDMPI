
#include "../gtest-1.6.0/include/gtest/gtest.h"
#include <mpi.h>
#include <cstdlib>

using ::testing::EmptyTestEventListener;
using ::testing::InitGoogleTest;
using ::testing::Test;
using ::testing::TestCase;
using ::testing::TestEventListeners;
using ::testing::TestInfo;
using ::testing::TestPartResult;
using ::testing::UnitTest;

enum GTestColor {
  COLOR_DEFAULT,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_YELLOW
};

class MPI_PrettyUnitTestResultPrinter : public ::testing::TestEventListener {
public:
  // Get basic MPI info.
  MPI_PrettyUnitTestResultPrinter() {
	  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	  MPI_Comm_size(MPI_COMM_WORLD, &np);
	  
  }

  static void PrintTestName(const char * test_case, const char * test) {
    printf("%s.%s", test_case, test);
  }

  // The following methods override what's in the TestEventListener class.
  virtual void OnTestProgramStart(const UnitTest& unit_test) {}
  virtual void OnTestIterationStart(const UnitTest& unit_test, int iteration);
  virtual void OnEnvironmentsSetUpStart(const UnitTest& unit_test);
  virtual void OnEnvironmentsSetUpEnd(const UnitTest& unit_test) {}
  virtual void OnTestCaseStart(const TestCase& test_case);
  virtual void OnTestStart(const TestInfo& test_info);
  virtual void OnTestPartResult(const TestPartResult& result);
  virtual void OnTestEnd(const TestInfo& test_info);
  virtual void OnTestCaseEnd(const TestCase& test_case);
  virtual void OnEnvironmentsTearDownStart(const UnitTest& unit_test);
  virtual void OnEnvironmentsTearDownEnd(const UnitTest& unit_test) {}
  virtual void OnTestIterationEnd(const UnitTest& unit_test, int iteration);
  virtual void OnTestProgramEnd(const UnitTest& unit_test) {}

private:
  int np, rank;

  static void PrintFailedTests(const UnitTest& unit_test);
  std::string test_case_name_;
};



