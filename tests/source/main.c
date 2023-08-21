/* ****************************************************************************************************************** */
/**
 *  @file main.c
 *  @author Edward J. Parkinson (e.parkinson@soton.ac.uk)
 *  @date August 2023
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

#include <stdlib.h>
#include <check.h>

#include "tests/tests.h"

#ifdef CUDA_ON
#include <cuda_runtime.h>
#endif

/** *******************************************************************************************************************
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

int
main (void)
{
  int number_failed;
  Suite *matrix_suite;
  SRunner *sr;

#ifdef CUDA_ON
  cudaSetDevice (0);
  cudaDeviceEnablePeerAccess (0, 0);
  cuda_init ();
#endif

  // Create the suite
  matrix_suite = create_matrix_suite ();

  // Create a test runner
  sr = srunner_create (matrix_suite);

  // Run the tests
  srunner_run_all (sr, CK_NORMAL);

  // Get the number of failed tests
  number_failed = srunner_ntests_failed (sr);

  // Clean up
  srunner_free (sr);

#ifdef CUDA_ON
  cuda_finish ();
#endif

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;    // Return 0 if all tests passed, 1 if there are failures
}
