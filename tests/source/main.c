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

int cuda_init (void);
int cuda_finish (void);



/** *******************************************************************************************************************
 *
 *  @brief
 *
 *  ***************************************************************************************************************** */

int
main (void)
{
#ifdef CUDA_ON
  cuda_init ();
#endif

  int number_failed;
  Suite *s;
  SRunner *sr;

  // Create the suite
  s = matrix_suite ();

  // Create a test runner
  sr = srunner_create (s);

  // Run the tests
  srunner_run_all (sr, CK_NORMAL);

  // Get the number of failed tests
  number_failed = srunner_ntests_failed (sr);

  // Clean up
  srunner_free (sr);

  return (number_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;    // Return 0 if all tests passed, 1 if there are failures

#ifdef CUDA_ON
  cuda_finish ();
#endif

  return EXIT_SUCCESS;
}
