/*
 * gtest_flow.cc
 *
 *  Created on: Mar 15, 2012
 *      Author: jb
 */

#include "gtest_mpi.h"

void MPI_PrettyUnitTestResultPrinter::OnTestIterationStart(const UnitTest& unit_test, int iteration) {
/*
    if (GTEST_FLAG(repeat) != 1)
    printf("\nRepeating all tests (iteration %d) . . .\n\n", iteration + 1);

  const char* const filter = GTEST_FLAG(filter).c_str();

  // Prints the filter if it's not *.  This reminds the user that some
  // tests may be skipped.
  if (! internal::String::CStringEquals(filter, kUniversalFilter)) {
    ColoredPrintf(COLOR_YELLOW,
                  "Note: %s filter = %s\n", GTEST_NAME_, filter);
  }

  if (internal::ShouldShard(kTestTotalShards, kTestShardIndex, false)) {
    const Int32 shard_index = Int32FromEnvOrDie(kTestShardIndex, -1);
    ColoredPrintf(COLOR_YELLOW,
                  "Note: This is test shard %d of %s.\n",
                  static_cast<int>(shard_index) + 1,
                  internal::posix::GetEnv(kTestTotalShards));
  }

  if (GTEST_FLAG(shuffle)) {
    ColoredPrintf(COLOR_YELLOW,
                  "Note: Randomizing tests' orders with a seed of %d .\n",
                  unit_test.random_seed());
  }

  ColoredPrintf(COLOR_GREEN,  "[==========] ");
  printf("Running %s from %s.\n",
         FormatTestCount(unit_test.test_to_run_count()).c_str(),
         FormatTestCaseCount(unit_test.test_case_to_run_count()).c_str());
  fflush(stdout);
*/
}


void MPI_PrettyUnitTestResultPrinter::OnEnvironmentsSetUpStart(const UnitTest& unit_test) {
  if (rank==0) {
      ColoredPrintf(COLOR_GREEN,  "[----------] ");
      printf("Global test environment set-up.\n");
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_PrettyUnitTestResultPrinter::OnTestCaseStart(const TestCase& test_case) {
  if (rank==0) {
      test_case_name_ = test_case.name();
      const std::string counts =
              FormatCountableNoun(test_case.test_to_run_count(), "test", "tests");

      ColoredPrintf(COLOR_GREEN, "[----------] ");
      printf("%s from %s", counts.c_str(), test_case_name_.c_str());
      if (test_case.type_param() == NULL) {
          printf("\n");
      } else {
          printf(", where TypeParam = %s\n", test_case.type_param());
      }
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_PrettyUnitTestResultPrinter::OnTestStart(const TestInfo& test_info) {
  if (rank == 0) {
      ColoredPrintf(COLOR_GREEN,  "[ RUN      ] ");
      PrintTestName(test_case_name_.c_str(), test_info.name());
      printf("\n");
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// Called after a failed assertion or a SUCCEED() invocation.
void MPI_PrettyUnitTestResultPrinter::OnTestPartResult(const TestPartResult& result) {
    // let print one process after other
    for(unsigned int i=0; i < np; i++) {
        if (rank == i && (result.type() != TestPartResult::kSuccess)) {
            printf("[%d] ", rank);
            PrintTestPartResult(result);
            fflush(stdout);
        }
      //  MPI_Barrier(MPI_COMM_WORLD);
    }
}

// (no yet fully MPI friendly)
void MPI_PrettyUnitTestResultPrinter::OnTestEnd(const TestInfo& test_info) {
  int success=(test_info.result()->Passed());
  int reduced_success;
  MPI_Allreduce(&success, &reduced_success, 1, MPI_INT, MPI_MIN,MPI_COMM_WORLD);

  if (rank == 0) {
      if (reduced_success) {
          ColoredPrintf(COLOR_GREEN, "[       OK ] ");
      } else {
          ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
      }
      PrintTestName(test_case_name_.c_str(), test_info.name());

      // HERE we should print comments for failed processes
      if (! reduced_success)
          PrintFullTestCommentIfPresent(test_info);

      //if (GTEST_FLAG(print_time)) {
      //    printf(" (%s ms)\n", internal::StreamableToString(
      //            test_info.result()->elapsed_time()).c_str());
      //} else {
          printf("\n");
	  //}
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_PrettyUnitTestResultPrinter::OnTestCaseEnd(const TestCase& test_case) {
  if (!GTEST_FLAG(print_time)) return;

  if (rank==0) {
      test_case_name_ = test_case.name();
      const std::string counts =
              FormatCountableNoun(test_case.test_to_run_count(), "test", "tests");
      ColoredPrintf(COLOR_GREEN, "[----------] ");
      //printf("%s from %s (%s ms total)\n\n",
      //   counts.c_str(), test_case_name_.c_str(),
      //   internal::StreamableToString(test_case.elapsed_time()).c_str());
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void MPI_PrettyUnitTestResultPrinter::OnEnvironmentsTearDownStart(
    const UnitTest& /*unit_test*/) {
  if (rank==0) {
      ColoredPrintf(COLOR_GREEN,  "[----------] ");
      printf("Global test environment tear-down\n");
      fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

// Internal helper for printing the list of failed tests.
// (no yet MPI friendly)
void MPI_PrettyUnitTestResultPrinter::PrintFailedTests(const UnitTest& unit_test) {
  const int failed_test_count = unit_test.failed_test_count();
  if (failed_test_count == 0) {
    return;
  }

  for (int i = 0; i < unit_test.total_test_case_count(); ++i) {
    const TestCase& test_case = *unit_test.GetTestCase(i);
    if (!test_case.should_run() || (test_case.failed_test_count() == 0)) {
      continue;
    }
    for (int j = 0; j < test_case.total_test_count(); ++j) {
      const TestInfo& test_info = *test_case.GetTestInfo(j);
      if (!test_info.should_run() || test_info.result()->Passed()) {
        continue;
      }
      ColoredPrintf(COLOR_RED, "[  FAILED  ] ");
      printf("%s.%s", test_case.name(), test_info.name());
      PrintFullTestCommentIfPresent(test_info);
      printf("\n");
    }
  }
}

// (no yet MPI friendly)
void MPI_PrettyUnitTestResultPrinter::OnTestIterationEnd(const UnitTest& unit_test,
							 int iteration) {
  ColoredPrintf(COLOR_GREEN,  "[==========] ");
  printf("%s from %s ran.",
         FormatTestCount(unit_test.test_to_run_count()).c_str(),
         FormatTestCaseCount(unit_test.test_case_to_run_count()).c_str());
  //if (GTEST_FLAG(print_time)) {
  //  printf(" (%s ms total)",
  //         internal::StreamableToString(unit_test.elapsed_time()).c_str());
  //}
  printf("\n");
  ColoredPrintf(COLOR_GREEN,  "[  PASSED  ] ");
  printf("%s.\n", FormatTestCount(unit_test.successful_test_count()).c_str());

  int num_failures = unit_test.failed_test_count();
  if (!unit_test.Passed()) {
    const int failed_test_count = unit_test.failed_test_count();
    ColoredPrintf(COLOR_RED,  "[  FAILED  ] ");
    printf("%s, listed below:\n", FormatTestCount(failed_test_count).c_str());
    PrintFailedTests(unit_test);
    printf("\n%2d FAILED %s\n", num_failures,
                        num_failures == 1 ? "TEST" : "TESTS");
  }

  int num_disabled = unit_test.disabled_test_count();
  if (num_disabled && !GTEST_FLAG(also_run_disabled_tests)) {
    if (!num_failures) {
      printf("\n");  // Add a spacer if no FAILURE banner is displayed.
    }
    ColoredPrintf(COLOR_YELLOW,
                  "  YOU HAVE %d DISABLED %s\n\n",
                  num_disabled,
                  num_disabled == 1 ? "TEST" : "TESTS");
  }
  // Ensure that Google Test output is printed before, e.g., heapchecker output.
  fflush(stdout);
}





